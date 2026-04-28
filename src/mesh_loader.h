#pragma once

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <unordered_map>
#include <filesystem>
#include <Eigen/Dense>

using namespace Eigen;

// =====================
// DATA STRUCTURES
// =====================

struct Material {
    std::string name;
    Vector3f kd = Vector3f(0.75f, 0.75f, 0.75f);
    std::string map_kd; // texture filename
};

struct TriMesh {
    std::vector<Vector3f> vertices;
    std::vector<Vector3i> faces;

    std::vector<Vector2f> texcoords;
    std::vector<Vector3i> face_texcoord_ids;

    std::vector<Material> materials;
    std::vector<int> face_material_ids;
};

// =====================
// HELPERS
// =====================

static int parseObjIndex(const std::string& token) {
    std::stringstream ss(token);
    std::string indexStr;
    std::getline(ss, indexStr, '/');
    return std::stoi(indexStr) - 1;
}

static int parseObjTexIndex(const std::string& token) {
    std::stringstream ss(token);
    std::string vStr, vtStr;

    std::getline(ss, vStr, '/');
    std::getline(ss, vtStr, '/');

    if (vtStr.empty()) return -1;
    return std::stoi(vtStr) - 1;
}

// =====================
// LOAD MTL
// =====================

inline std::vector<Material> load_mtl(
    const std::string& mtl_path,
    std::unordered_map<std::string, int>& mat_name_to_id
    ) {
    std::vector<Material> materials;

    std::ifstream file(mtl_path);
    if (!file.is_open()) {
        std::cout << "Could not open mtl: " << mtl_path << std::endl;
        return materials;
    }

    std::string line;
    Material current;
    bool has_material = false;

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string type;

        if (!(ss >> type)) continue;

        if (type == "newmtl") {
            if (has_material) {
                mat_name_to_id[current.name] = (int)materials.size();
                materials.push_back(current);
            }

            current = Material();
            ss >> current.name;
            has_material = true;
        }
        else if (type == "Kd") {
            ss >> current.kd.x() >> current.kd.y() >> current.kd.z();
        }
        else if (type == "map_Kd") {
            ss >> current.map_kd;
        }
    }

    if (has_material) {
        mat_name_to_id[current.name] = (int)materials.size();
        materials.push_back(current);
    }

    return materials;
}

// =====================
// LOAD OBJ
// =====================

inline TriMesh load_obj(const std::string& filepath) {
    TriMesh mesh;

    std::ifstream file(filepath);
    if (!file.is_open()) {
        std::cout << "Could not open obj: " << filepath << std::endl;
        return mesh;
    }

    std::unordered_map<std::string, int> mat_name_to_id;
    int current_material_id = -1;

    std::filesystem::path obj_path(filepath);
    std::string obj_dir = obj_path.parent_path().string();

    std::string line;

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string type;

        if (!(ss >> type)) continue;

        // -----------------
        // VERTICES
        // -----------------
        if (type == "v") {
            float x, y, z;
            ss >> x >> y >> z;
            mesh.vertices.push_back(Vector3f(x, y, z));
        }

        // -----------------
        // UVs
        // -----------------
        else if (type == "vt") {
            float u, v;
            ss >> u >> v;
            mesh.texcoords.push_back(Vector2f(u, v));
        }

        // -----------------
        // MTL FILE
        // -----------------
        else if (type == "mtllib") {
            std::string mtl_file;
            ss >> mtl_file;

            std::string mtl_path =
                obj_dir.empty() ? mtl_file : obj_dir + "/" + mtl_file;

            mesh.materials = load_mtl(mtl_path, mat_name_to_id);
        }

        // -----------------
        // MATERIAL SWITCH
        // -----------------
        else if (type == "usemtl") {
            std::string mat_name;
            ss >> mat_name;

            if (mat_name_to_id.count(mat_name)) {
                current_material_id = mat_name_to_id[mat_name];
            } else {
                current_material_id = -1;
                std::cout << "Unknown material: " << mat_name << std::endl;
            }
        }

        // -----------------
        // FACES
        // -----------------
        else if (type == "f") {
            std::vector<std::string> tokens;
            std::string tok;

            while (ss >> tok) {
                tokens.push_back(tok);
            }

            if (tokens.size() < 3) {
                continue;
            }

            // Fan triangulation:
            // quad: 0,1,2 and 0,2,3
            // ngon: 0,1,2 then 0,2,3 then 0,3,4 ...
            for (int t = 1; t + 1 < tokens.size(); t++) {
                int a = parseObjIndex(tokens[0]);
                int b = parseObjIndex(tokens[t]);
                int c = parseObjIndex(tokens[t + 1]);

                int ta = parseObjTexIndex(tokens[0]);
                int tb = parseObjTexIndex(tokens[t]);
                int tc = parseObjTexIndex(tokens[t + 1]);

                mesh.faces.push_back(Vector3i(a, b, c));
                mesh.face_texcoord_ids.push_back(Vector3i(ta, tb, tc));
                mesh.face_material_ids.push_back(current_material_id);
            }
        }
    }

    std::cout << "Loaded OBJ: " << filepath << std::endl;
    std::cout << "verts: " << mesh.vertices.size()
              << " faces: " << mesh.faces.size()
              << " materials: " << mesh.materials.size()
              << " texcoords: " << mesh.texcoords.size()
              << std::endl;

    return mesh;
}
