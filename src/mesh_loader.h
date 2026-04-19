#pragma once

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <cassert>
#include <iostream>
#include <Eigen/Dense>

using namespace Eigen;

struct TriMesh {
    std::vector<Vector3f> vertices;
    std::vector<Vector3i> faces;
};


inline TriMesh load_obj(const std::string& filepath, float scale = 1.0f, bool center = false) {
    std::vector<Vector3f> verts;
    std::vector<Vector3i> faces;

    std::ifstream fh(filepath);
    if (!fh.is_open()) {
        std::cerr << "mesh_loader: failed to open " << filepath << std::endl;
        return {};
    }

    std::string line;
    while (std::getline(fh, line)) {
        if (line.empty() || line[0] == '#')
            continue;

        std::istringstream ss(line);
        std::string key;
        ss >> key;

        if (key == "v") {
            float x, y, z;
            ss >> x >> y >> z;
            verts.push_back(Vector3f(x, y, z));

        } else if (key == "f") {
            std::vector<int> face_vi;
            std::string token;
            while (ss >> token) {
                int vi = std::stoi(token.substr(0, token.find('/'))) - 1;
                face_vi.push_back(vi);
            }

            // Fan triangulate for quads and ngons
            for (int k = 1; k < (int)face_vi.size() - 1; k++) {
                faces.push_back(Vector3i(face_vi[0], face_vi[k], face_vi[k + 1]));
            }
        }
        // vn and vt lines skipped — normals recomputed from geometry
    }

    // Apply scale (e.g. 0.01f to convert cm -> m)
    for (auto& v : verts) {
        v *= scale;
    }

    // Center at origin if requested
    if (center && !verts.empty()) {
        Vector3f centroid = Vector3f::Zero();
        for (const auto& v : verts) {
            centroid += v;
        }
        centroid /= (float)verts.size();
        for (auto& v : verts) {
            v -= centroid;
        }
    }

    TriMesh mesh;
    mesh.vertices = verts;
    mesh.faces = faces;
    return mesh;
}


inline std::vector<TriMesh> load_pose_sequence(const std::vector<std::string>& filepaths, float scale = 1.0f, bool center = false) {
    std::vector<TriMesh> meshes;
    for (const auto& fp : filepaths) {
        meshes.push_back(load_obj(fp, scale, center));
    }

    assert(!meshes.empty());

    const auto& ref = meshes[0];
    for (int i = 1; i < (int)meshes.size(); i++) {
        assert(meshes[i].vertices.size() == ref.vertices.size());
        assert(meshes[i].faces.size() == ref.faces.size());
    }

    return meshes;
}
