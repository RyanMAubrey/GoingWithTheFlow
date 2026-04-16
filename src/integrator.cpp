#include "integrator.h"
#include "momentum.h"
#include "lift_and_drag.cpp"
#include "mesh_loader.h"

#include <set>

Integrator::Integrator() {}

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

    // Save unique edges into final edge vector
    for (std::pair<int,int> e : unique_edges) {
        // TODO: Figure out how to calculate bending angles
        edges.push_back({e.first, e.second, 0.0f});
    }
}

void Integrator::CalculateFaceAttributes(std::vector<Vector3i>& faces, std::vector<Vector3f>& gamma) {
    face_areas.resize(faces.size());
    face_normals.resize(faces.size());

    for (int f = 0; f < faces.size(); f++) {
        int i = faces[f].x();
        int j = faces[f].y();
        int k = faces[f].z();
        Vector3f normal = gamma[j] - gamma[i].cross(gamma[k] - gamma[i]);

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
