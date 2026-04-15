#include "momentum.h"

using namespace Eigen;
using Vector6f = Eigen::Matrix<float, 6, 1>;

Vector6f Momentum::make_vec6(const Vector3f& rot, const Vector3f& trans) {
    Vector6f v;
    v << rot, trans;
    return v;
}

Momentum::Momentum() {}

// h = time step
// gamma_k = all vertex positions at time k
// gamma_k1 = all vertex positions at time k+1
// mass density = non-negative
Vector6f Momentum::calc_body_momentum(std::vector<Vector3f>& gamma_k, std::vector<Vector3f>& gamma_k1,
                                       std::vector<float>& mass_density, float h) {
    Vector3f l_b(0.0f);
    Vector3f p_b(0.0f);

    for (int i = 0; i < gamma_k.size(); i++) {
        float rho_b = mass_density[i];
        Vector3f gamma_prime_i = (gamma_k1[i] - gamma_k[i]) / h;

        l_b += gamma_k[i].cross(gamma_prime_i) * rho_b;
        p_b += gamma_prime_i * rho_b;
    }

    return make_vec6(l_b, p_b);
}

// h = time step
// gamma_k = all vertex positions at time k
// gamma_k1 = all vertex positions at time k+1
// faces = each face has i,j,k
// edges = edge struct
Vector6f Momentum::calc_fluid_momentum(std::vector<Vector3f>& gamma_k, std::vector<Vector3f>& gamma_k1,
                         std::vector<Vector3i>& faces, std::vector<Edge>& edges,
                         float rho_f, float h) {
    // Calculate face values
    std::vector<float> face_areas(faces.size());
    std::vector<Vector3f> face_normals(faces.size());
    for (int f = 0; f < faces.size(); f++) {
        int i = faces[f].x();
        int j = faces[f].y();
        int k = faces[f].z();
        Vector3f normal = gamma_k[j] - gamma_k[i].cross(gamma_k[k] - gamma_k[i]);

        face_areas[f] = 0.5f * normal.norm(); // area = 1/2 * ||N||
        face_normals[f] = normal.normalized();
    }

    // Calculate delta
    float total_area = 0.0f;
    for (int f = 0; f < faces.size(); f++) {
        total_area += face_areas[f];
    }
    float total_edge = 0.0f;
    for (int e = 0; e < edges.size(); e++) {
        // Given edge index, find vertices for length & bending angle
        Edge curr_e = edges[e];
        float alpha_ij = curr_e.alpha;
        float l_ij = (gamma_k[curr_e.j] - gamma_k[curr_e.i]).norm();
        total_edge += alpha_ij * l_ij;
    }
    float delta = total_area / total_edge;

    Vector3f l_f(0.0f);
    Vector3f p_f(0.0f);
    for (int f = 0; f < faces.size(); f++) {
        int i = faces[f].x();
        int j = faces[f].y();
        int k = faces[f].z();
        // Combine ijk into one value (centers)
        Vector3f gamma_ijk = (gamma_k[i] + gamma_k[j] + gamma_k[k]) / 3.0f;
        Vector3f gamma_prime_ijk = ((gamma_k1[i] - gamma_k[i]) + (gamma_k1[j] - gamma_k[j]) + (gamma_k1[k] - gamma_k[k])) / (3.0f * h);

        Vector3f n_ijk = face_normals[f];
        float A_ijk = face_areas[f];
        float dot_prod = gamma_prime_ijk.dot(n_ijk);

        l_f += rho_f * delta * dot_prod * gamma_ijk.cross(n_ijk) * A_ijk;
        p_f += rho_f * delta * dot_prod * n_ijk * A_ijk;
    }

    return make_vec6(l_f, p_f);
}
