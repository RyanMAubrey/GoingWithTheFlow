#include "integrator.h"
#include "momentum.h"
#include "inertia.h"
#include "se3.h"

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
                          const Vector3f& background_flow,
                          float drag_scale);

float compute_volume(const std::vector<Vector3f>& gamma,
                     const std::vector<Vector3i>& faces);

#include <iostream>
#include <set>
#include <algorithm>
#include <fstream>

Integrator::Integrator() {}

void Integrator::Simulate() {
    Simulate(AddedMassMode::Global);
}

void Integrator::Simulate(AddedMassMode mode) {
    const char* mode_name = (mode == AddedMassMode::Piecewise) ? "piecewise" : "global";
    std::cout << "[simulate:" << mode_name << "] starting" << std::endl;
    output_frames.clear();
    output_centroids.clear();
    LoadAllPoses();

    // Initial state
    Matrix4f pose = Matrix4f::Identity();
    Vector6f body_velocity = Vector6f::Zero();
    Vector6f mu = Vector6f::Zero();
    
    // Compute per-vertex mass from density * volume / num_vertices
    float volume = compute_volume(all_vertices[0], shared_faces);
    float total_body_mass = rho_f * volume;
    float per_vertex_mass = total_body_mass / (float)all_vertices[0].size();
    std::vector<float> mass_density(all_vertices[0].size(), per_vertex_mass);

    Momentum m;

    int global_step = 0;
    int completed_pairs = 0;
    const int total_pairs = num_strokes * (total_frames - 1);
    const int progress_interval = std::max(1, total_pairs / 10);
    bool printed_delta_summary = false;

    for (int swim_stroke = 0; swim_stroke < num_strokes; swim_stroke++) {
        std::cout << "[simulate:" << mode_name << "] stroke "
                  << (swim_stroke + 1) << "/" << num_strokes << std::endl;
        for (int i = 0; i < total_frames - 1; i++) {

            std::vector<Vector3f>& frame_A = all_vertices[i];
            std::vector<Vector3f>& frame_B = all_vertices[i + 1];
            int num_verts = (int)frame_A.size();

            for (int substep = 0; substep < substeps; substep++) {
                float alpha = (float) substep / (float) substeps;
                float alpha_next = (float) (substep + 1) / (float) substeps;

                // Interpolate positions
                std::vector<Vector3f> vertices_k(num_verts), vertices_k1(num_verts);
                for (int v = 0; v < num_verts; v++) {
                    vertices_k[v]  = (1.0f - alpha) * frame_A[v] + alpha * frame_B[v];
                    vertices_k1[v] = (1.0f - alpha_next) * frame_A[v] + alpha_next * frame_B[v];
                }

            // Recompute face attributes for this frame
            std::vector<float> face_areas;
            std::vector<Vector3f> face_normals;
            CalculateFaceAttributes(shared_faces, vertices_k, face_areas, face_normals);
            
            float delta = CalculateDelta(shared_faces, shared_edges, vertices_k, face_areas);
            std::vector<float> face_deltas;
            std::array<float, NUM_REGIONS> region_deltas = {};
            if (mode == AddedMassMode::Piecewise &&
                shared_face_regions.size() == shared_faces.size()) {
                region_deltas = CalculateRegionDeltas(vertices_k, face_areas, delta);
                face_deltas = BuildFaceDeltas(region_deltas, delta);
            }

            if (!printed_delta_summary) {
                std::cout << "[simulate:" << mode_name << "] global delta: " << delta << std::endl;
                if (mode == AddedMassMode::Piecewise && !face_deltas.empty()) {
                    static const char* region_names[NUM_REGIONS] = {
                        "head",
                        "left_front_flipper",
                        "right_front_flipper",
                        "left_rear_flipper",
                        "right_rear_flipper",
                        "shell"
                    };
                    for (int region = 0; region < NUM_REGIONS; region++) {
                        std::cout << "[simulate:" << mode_name << "] "
                                  << region_names[region] << " delta: "
                                  << region_deltas[region] << std::endl;
                    }
                }
                printed_delta_summary = true;
            }
            
            // ————— Algorithm 1 —————
            // Lines 1-6: Kirchhoff Tensor and Momentum
            Matrix6f body_inertia  = calc_body_inertia(vertices_k, mass_density);
            Matrix6f added_mass =
                (mode == AddedMassMode::Piecewise && !face_deltas.empty())
                    ? calc_added_mass_piecewise(vertices_k, shared_faces, face_areas, face_normals, rho_f, face_deltas)
                    : calc_added_mass(vertices_k, shared_faces, face_areas, face_normals, rho_f, delta);
            Matrix6f kirchhoff_tensor = body_inertia + added_mass;

            Vector6f body_momentum = m.calc_body_momentum(vertices_k, vertices_k1, mass_density, h);
            Vector6f fluid_momentum =
                (mode == AddedMassMode::Piecewise && !face_deltas.empty())
                    ? m.calc_fluid_momentum_piecewise(vertices_k, vertices_k1, shared_faces, face_areas, face_normals, rho_f, h, face_deltas)
                    : m.calc_fluid_momentum(vertices_k, vertices_k1, shared_faces, face_areas, face_normals, shared_edges, rho_f, h, delta);
            Vector6f combined_momentum = body_momentum + fluid_momentum;

            Vector6f total_forces = calc_total_force(vertices_k, vertices_k1, shared_faces,
                                          pose, body_velocity, rho_f, h, gravity, total_mass, bg_flow, drag_scale);
            
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

            // Added to prevent turtle from falling back
            body_velocity.head<3>() *= angular_damping;

            mu = kirchhoff_tensor * body_velocity + combined_momentum;
            pose = pose * cayley_map(h * body_velocity);

            // Position mesh in world space and save
            Matrix3f A = pose.block<3,3>(0,0);
            Vector3f b = pose.block<3,1>(0,3);

            std::vector<Vector3f> world_verts(vertices_k1.size());
            for (int v = 0; v < (int)vertices_k1.size(); v++) {
                world_verts[v] = A * vertices_k1[v] + b;
            }

            // Store current frame for visualizer
            if (global_step % output_every == 0) {
                std::vector<Eigen::Vector3d> frame_d(world_verts.size());
                Eigen::Vector3d centroid = Eigen::Vector3d::Zero();
                for (int v = 0; v < (int)world_verts.size(); v++) {
                    frame_d[v] = world_verts[v].cast<double>();
                    centroid += frame_d[v];
                }
                centroid /= std::max(1.0, static_cast<double>(world_verts.size()));
                output_frames.push_back(std::move(frame_d));
                output_centroids.push_back(centroid);
            }

            global_step++;
            }

            completed_pairs++;
            if (completed_pairs % progress_interval == 0 || completed_pairs == total_pairs) {
                std::cout << "[simulate:" << mode_name << "] "
                          << completed_pairs << "/" << total_pairs
                          << " pose pairs complete" << std::endl;
            }
        }
    }

    std::cout << "[simulate:" << mode_name << "] complete, stored "
              << output_frames.size() << " output frames" << std::endl;
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

bool Integrator::HasEdge(const Vector3i& face, int i, int j) const {
    bool has_i = (face[0] == i || face[1] == i || face[2] == i);
    bool has_j = (face[0] == j || face[1] == j || face[2] == j);
    return has_i && has_j;
}

void Integrator::BuildEdgeFaceAdjacency() {
    shared_edge_faces.clear();
    shared_edge_faces.reserve(shared_edges.size());

    for (const Edge& edge : shared_edges) {
        int face_a = -1;
        int face_b = -1;

        for (int f = 0; f < (int)shared_faces.size(); f++) {
            if (!HasEdge(shared_faces[f], edge.i, edge.j)) {
                continue;
            }

            if (face_a == -1) {
                face_a = f;
            } else {
                face_b = f;
                break;
            }
        }

        shared_edge_faces.push_back({face_a, face_b});
    }
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

std::array<float, NUM_REGIONS> Integrator::CalculateRegionDeltas(const std::vector<Vector3f>& gamma,
                                                                 const std::vector<float>& face_areas,
                                                                 float fallback_delta) const {
    std::array<float, NUM_REGIONS> region_areas = {};
    std::array<float, NUM_REGIONS> region_edge_sums = {};
    std::array<float, NUM_REGIONS> region_deltas = {};

    for (int f = 0; f < (int)shared_face_regions.size() && f < (int)face_areas.size(); f++) {
        const int region = shared_face_regions[f];
        if (region >= 0 && region < NUM_REGIONS) {
            region_areas[region] += face_areas[f];
        }
    }

    for (int e = 0; e < (int)shared_edges.size() && e < (int)shared_edge_faces.size(); e++) {
        const Edge& edge = shared_edges[e];
        const std::pair<int, int>& edge_faces = shared_edge_faces[e];
        const float contribution = edge.alpha * (gamma[edge.j] - gamma[edge.i]).norm();

        const int face_a = edge_faces.first;
        const int face_b = edge_faces.second;
        if (face_a >= 0 && face_a < (int)shared_face_regions.size()) {
            const int region_a = shared_face_regions[face_a];
            if (region_a >= 0 && region_a < NUM_REGIONS) {
                region_edge_sums[region_a] += contribution;
            }
        }
        if (face_b >= 0 && face_b < (int)shared_face_regions.size()) {
            const int region_b = shared_face_regions[face_b];
            if (region_b >= 0 && region_b < NUM_REGIONS &&
                (face_a < 0 || region_b != shared_face_regions[face_a])) {
                region_edge_sums[region_b] += contribution;
            }
        }
    }

    for (int region = 0; region < NUM_REGIONS; region++) {
        if (region_areas[region] > 1e-8f && region_edge_sums[region] > 1e-8f) {
            region_deltas[region] = region_areas[region] / region_edge_sums[region];
        } else {
            region_deltas[region] = fallback_delta;
        }
    }

    return region_deltas;
}

std::vector<float> Integrator::BuildFaceDeltas(const std::array<float, NUM_REGIONS>& region_deltas,
                                               float fallback_delta) const {
    std::vector<float> face_deltas(shared_faces.size(), fallback_delta);
    for (int f = 0; f < (int)shared_faces.size() && f < (int)shared_face_regions.size(); f++) {
        const int region = shared_face_regions[f];
        if (region >= 0 && region < NUM_REGIONS) {
            face_deltas[f] = region_deltas[region];
        }
    }
    return face_deltas;
}

// Pose Loading
void Integrator::LoadAllPoses() {
    all_vertices.clear();
    all_vertices.resize(total_frames);
    shared_face_regions.clear();
    shared_edge_faces.clear();

    for (int i = 0; i < total_frames; i++) {
        std::string path = animal + "_poses/frame_" + std::to_string(i) + ".obj";

        if (i == 0) {
            // First frame: load full topology (faces, edges).
            TriMesh mesh = LoadPose(path, all_vertices[i], shared_faces, shared_edges,
                     shared_face_areas_scratch, shared_face_normals_scratch);

            // Store texture data for rendering
            texcoords = mesh.texcoords;
            face_texcoord_ids = mesh.face_texcoord_ids;
            if (!mesh.materials.empty() && !mesh.materials[0].map_kd.empty()) {
                texture_path = animal + "_poses/" + mesh.materials[0].map_kd;
            }
            BuildEdgeFaceAdjacency();
            shared_face_regions = partition_mesh(all_vertices[i], shared_faces);
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
    }
    std::cout << "All " << total_frames << " poses loaded." << std::endl;
}
