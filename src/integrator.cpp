#include "integrator.h"
#include "momentum.h"
#include "inertia.h"
// Forward declarations from lift_and_drag.cpp
Vector6f calc_total_force(const std::vector<Vector3f>& positions_k,
                          const std::vector<Vector3f>& positions_k1,
                          const std::vector<Vector3i>& faces,
                          const Matrix4f& pose,
                          const Vector6f& body_twist,
                          float fluid_density,
                          float timestep,
                          const Vector3f& gravity,
                          float total_mass,
                          const Vector3f& background_flow);
#include "se3.h"

#include <iostream>
#include <set>
#include <algorithm>
#include <fstream>

Integrator::Integrator() {}

static void write_obj(const std::string& path,
                      const std::vector<Vector3f>& verts,
                      const std::vector<Vector3i>& faces) {
    std::ofstream file(path);
    if (!file.is_open()) {
        std::cerr << "Could not write: " << path << std::endl;
        return;
    }
    for (const auto& v : verts) {
        file << "v " << v.x() << " " << v.y() << " " << v.z() << "\n";
    }
    for (const auto& f : faces) {
        file << "f " << (f.x()+1) << " " << (f.y()+1) << " " << (f.z()+1) << "\n";
    }
}

void Integrator::Simulate() {
    LoadAllPoses();

    // Initial state
    Matrix4f pose = Matrix4f::Identity();
    Vector6f body_velocity = Vector6f::Zero();
    Vector6f mu = Vector6f::Zero();

    // Scene parameters
    float total_mass = 0.0f;                    // neutrally buoyant
    Vector3f gravity(0.0f, -9.81f, 0.0f);       // standard gravity
    Vector3f bg_flow = Vector3f::Zero();        // no background flow for now
    std::vector<float> mass_density(all_vertices[0].size(), rho_f);  // same as water

    Momentum m;

    for (int swim_stroke = 0; swim_stroke < num_strokes; swim_stroke++) {
        for (int i = 0; i < total_frames - 1; i++) {
            int step = swim_stroke * (total_frames - 1) + i;
            
            std::vector<Vector3f>& vertices_k  = all_vertices[i];
            std::vector<Vector3f>& vertices_k1 = all_vertices[i + 1];

            // Recompute face attributes for this frame
            std::vector<float> face_areas;
            std::vector<Vector3f> face_normals;
            CalculateFaceAttributes(shared_faces, vertices_k, face_areas, face_normals);
            
            float delta = CalculateDelta(shared_faces, shared_edges, vertices_k, face_areas);
            
            // ————— Algorithm 1 —————
            // Lines 1-6: Kirchhoff Tensor and Momentum
            Matrix6f body_inertia  = calc_body_inertia(vertices_k, mass_density);
            Matrix6f added_mass = calc_added_mass(vertices_k, shared_faces, face_areas, face_normals, rho_f, delta);
            Matrix6f kirchhoff_tensor = body_inertia + added_mass;

            Vector6f body_momentum = m.calc_body_momentum(vertices_k, vertices_k1, mass_density, h);
            Vector6f fluid_momentum = m.calc_fluid_momentum(vertices_k, vertices_k1, shared_faces, face_areas, face_normals, shared_edges, rho_f, h, delta);
            Vector6f combined_momentum = body_momentum + fluid_momentum;

            Vector6f total_forces = calc_total_force(vertices_k, vertices_k1, shared_faces,
                                          pose, body_velocity, rho_f, h, gravity, total_mass, bg_flow);
            
            // Line 7: Newton Solve
            Vector6f rhs = dtau_inv_star(-h * body_velocity) * mu + h * total_forces;   // RHS is constant

            Vector6f body_velocity_new = body_velocity;  // initial guess
            for (int iter = 0; iter < 10; iter++) {
                Vector6f lhs = dtau_inv_star(-h * body_velocity_new) * 
                                (kirchhoff_tensor * body_velocity_new + combined_momentum);
                Vector6f residual = lhs - rhs;

                if (residual.norm() < 1e-8f) break;

                // Finite-difference Jacobian
                Matrix6f J;
                float eps = 1e-5f;
                for (int col = 0; col < 6; col++) {
                    Vector6f body_velocity_pert = body_velocity_new;
                    body_velocity_pert[col] += eps;
                    Vector6f lhs_pert = dtau_inv_star(-h * body_velocity_pert) * 
                                        (kirchhoff_tensor * body_velocity_pert + combined_momentum);
                    J.col(col) = (lhs_pert - lhs) / eps;
                }

                body_velocity_new -= J.fullPivLu().solve(residual);
            }

            // Lines 8-10: Update State
            body_velocity  = body_velocity_new;
            mu = kirchhoff_tensor * body_velocity + combined_momentum;
            pose = pose * cayley_map(h * body_velocity);

            // Position mesh in world space and save
            Matrix3f A = pose.block<3,3>(0,0);
            Vector3f b = pose.block<3,1>(0,3);

            std::vector<Vector3f> world_verts(vertices_k1.size());
            for (int v = 0; v < (int)vertices_k1.size(); v++) {
                world_verts[v] = A * vertices_k1[v] + b;
            }

            // Write output OBJ
            std::string out_path = "output/frame_" + std::to_string(step) + ".obj";
            write_obj(out_path, world_verts, shared_faces);

            // Print progress
            std::cout << "Step " << step
            << " stroke = " << swim_stroke
            << " pos = (" << b.x() << ", " << b.y() << ", " << b.z() << ")"
            << " |v| = " << body_velocity.tail<3>().norm()
            << std::endl;
        }
    }
}

TriMesh Integrator::LoadPose(const std::string& filepath, std::vector<Vector3f>& vertices, std::vector<Vector3i>& faces, std::vector<Edge>& edges,
                            std::vector<float>& face_areas, std::vector<Vector3f>& face_normals) {
    TriMesh m0 = load_obj(filepath);
    vertices = m0.vertices;
    faces = m0.faces;

    // Get unique edges using a set (auto removes duplicate)
    std::set<std::pair<int,int>> unique_edges;
    for (const Vector3i& face : faces) {
        int a = face[0];
        int b = face[1];
        int c = face[2];

        // Sort edges (make them all have smaller index first)
        unique_edges.insert({std::min(a,b), std::max(a,b)});
        unique_edges.insert({std::min(b,c), std::max(b,c)});
        unique_edges.insert({std::min(a,c), std::max(a,c)});
    }

    // Compute face attributes
    CalculateFaceAttributes(faces, vertices, face_areas, face_normals);

    // Save unique edges into final edge vector
    for (std::pair<int,int> e : unique_edges) {
        int i = e.first;
        int j = e.second;

        // For bending angles, need to be adjacent triangles (share i & j)
        int f1 = -1;
        int f2 = -1;
        for (int f = 0; f < faces.size(); f++) {
            if (HasEdge(faces[f], i, j)) {
                // Need to find two triangles (move on once both are found)
                if (f1 == -1) {
                    f1 = f;
                } else {
                    f2 = f;
                    break;
                }
            }
        }
        float alpha = 0.0f;
        if (f1 != -1 && f2 != -1) {
            // Dot product since angles depend on direction of normals
            float d = face_normals[f1].dot(face_normals[f2]);
            alpha = std::acos(std::max(-1.0f, std::min(1.0f, d))); // Clamp to make sure acos works
        }
        edges.push_back({i, j, alpha});
    }
    return m0;
}

bool Integrator::HasEdge(Vector3i& face, int i, int j) {
    bool has_i = (face[0] == i || face[1] == i || face[2] == i);
    bool has_j = (face[0] == j || face[1] == j || face[2] == j);
    return has_i && has_j;
}

void Integrator::CalculateFaceAttributes(std::vector<Vector3i>& faces, std::vector<Vector3f>& gamma,
                                        std::vector<float>& face_areas, std::vector<Vector3f>& face_normals) {
    face_areas.resize(faces.size());
    face_normals.resize(faces.size());

    for (int f = 0; f < faces.size(); f++) {
        int i = faces[f].x();
        int j = faces[f].y();
        int k = faces[f].z();
        Vector3f normal = (gamma[j] - gamma[i]).cross(gamma[k] - gamma[i]);

        face_areas[f] = 0.5f * normal.norm(); // area = 1/2 * ||N||
        face_normals[f] = normal.normalized();
    }
}

float Integrator::CalculateDelta(std::vector<Vector3i>& faces, std::vector<Edge>& edges, std::vector<Vector3f>& gamma, std::vector<float>& face_areas) {
    float total_area = 0.0f;
    for (int f = 0; f < faces.size(); f++) {
        total_area += face_areas[f];
    }
    float total_edge = 0.0f;
    for (int e = 0; e < edges.size(); e++) {
        // Given edge index, find vertices for length & bending angle
        Edge curr_e = edges[e];
        float alpha_ij = curr_e.alpha;
        float l_ij = (gamma[curr_e.j] - gamma[curr_e.i]).norm();
        total_edge += alpha_ij * l_ij;
    }
    return total_area / total_edge;
}

// Pose Loading
void Integrator::LoadAllPoses() {
    all_vertices.clear();
    all_vertices.resize(total_frames);

    for (int i = 0; i < total_frames; i++) {
        int frame_number = i * 8;  
        std::string path = "turtle_poses/frame_" + std::to_string(frame_number) + ".obj";

        if (i == 0) {
            // First frame: load full topology (faces, edges).
            LoadPose(path, all_vertices[i], shared_faces, shared_edges,
                     shared_face_areas_scratch, shared_face_normals_scratch);
        } else {
            // Subsequent frames: only need vertex positions.
            std::vector<Vector3i> unused_faces;
            std::vector<Edge> unused_edges;
            std::vector<float> unused_areas;
            std::vector<Vector3f> unused_normals;
            LoadPose(path, all_vertices[i], unused_faces, unused_edges,
                     unused_areas, unused_normals);
        }

        // Convert centimeters to meters
        for (auto& v : all_vertices[i]) {
            v *= 0.01f;
        }

        std::cout << "Loaded frame " << frame_number
                  << ": " << all_vertices[i].size() << " vertices." << std::endl;
    }

    std::cout << "All " << total_frames << " poses loaded." << std::endl;
}
