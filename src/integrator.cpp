#include "integrator.h"
#include "momentum.h"
#include "lift_and_drag.cpp"
#include "mesh_loader.h"

#include <set>

Integrator::Integrator() {}

void Integrator::Simulate() {
    // TODO: Algo 1 here
}

void Integrator::LoadPose(const std::string& filepath, std::vector<Vector3f>& vertices, std::vector<Vector3i>& faces, std::vector<Edge>& edges) {
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
    CalculateFaceAttributes(faces, vertices);

    // Save unique edges into final edge vector
    for (std::pair<int,int> e : unique_edges) {
        int i = e.first;
        int j = e.second;

        // For bending angles, need to be adjacent triangles (share i & j)
        int f1, f2 = -1;
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

void Integrator::CalculateFaceAttributes(std::vector<Vector3i>& faces, std::vector<Vector3f>& gamma) {
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

float Integrator::CalculateDelta(std::vector<Vector3i>& faces, std::vector<Edge>& edges, std::vector<Vector3f>& gamma) {
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
