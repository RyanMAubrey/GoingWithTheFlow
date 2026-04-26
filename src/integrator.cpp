#include "integrator.h"
#include "momentum.h"
#include "inertia.h"
//#include "lift_and_drag.cpp"
#include "mesh_loader.h"

#include <iostream>
#include <set>
#include <algorithm>

Integrator::Integrator() {}

void Integrator::Simulate() {
    // TODO: Algo 1 here
    Momentum m;

    for (int i = 1; i < total_poses; i++) {
        // Variable set up for k and k+1
        std::vector<Vector3f> vertices_k, vertices_k1;
        std::vector<Vector3i> faces_k, faces_k1;
        std::vector<Edge> edges_k, edges_k1;
        std::vector<float> face_areas_k, face_areas_k1;
        std::vector<Vector3f> face_normals_k, face_normals_k1;
        std::string filepath_k = "turtle_poses/frame_" + std::to_string(i) + ".obj";
        std::string filepath_k1 = "turtle_poses/frame_" + std::to_string(i+1) + ".obj";

        // Load the two consecutive poses
        LoadPose(filepath_k, vertices_k, faces_k, edges_k, face_areas_k, face_normals_k);
        LoadPose(filepath_k1, vertices_k1, faces_k1, edges_k1, face_areas_k1, face_normals_k1);

        // Extra variable calcs
        float delta = CalculateDelta(faces_k, edges_k, vertices_k, face_areas_k);
        std::vector<float> mass_density(vertices_k.size(), 1.0f); // Start with uniform mass (may change later)

        Matrix6f body_inertia = calc_body_inertia(vertices_k, mass_density);
        Matrix6f added_mass = calc_added_mass(vertices_k, faces_k,face_areas_k, face_normals_k, rho_f, delta);
        Matrix6f kirchhoff_tensor = body_inertia + added_mass;

        Vector6f body_momentum = m.calc_body_momentum(vertices_k, vertices_k1, mass_density, h);
        Vector6f fluid_momentum = m.calc_fluid_momentum(vertices_k, vertices_k1, faces_k, face_areas_k, face_normals_k, edges_k, rho_f, h, delta);
        Vector6f combined_momentum = body_momentum + fluid_momentum;

    }
}

void Integrator::LoadPose(const std::string& filepath, std::vector<Vector3f>& vertices, std::vector<Vector3i>& faces, std::vector<Edge>& edges,
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
