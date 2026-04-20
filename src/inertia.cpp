#include "inertia.h"

// Skew-symmetric (hat) matrix: hat(a) * b == a.cross(b)
static Matrix3f hat(const Vector3f& a) {
    Matrix3f m;
    m <<    0, -a.z(),  a.y(),
         a.z(),     0, -a.x(),
        -a.y(),  a.x(),     0;
    return m;
}

// Algorithm 3: Body inertia matrix M_b
//
// Derived from calc_body_momentum: for rigid motion gamma_prime_i = omega x gamma_i + v,
// substituting gives M_b * [omega; v] = [l_b; p_b] where:
//
//   top-left  (3x3): I = sum_i rho[i] * (|g|^2 * I3 - g * g^T)   (inertia tensor)
//   top-right (3x3): hat(p_first)    where p_first = sum_i rho[i] * gamma_i
//   bot-left  (3x3): hat(p_first)^T  = -hat(p_first)
//   bot-right (3x3): m * I3
Matrix6f calc_body_inertia(const std::vector<Vector3f>& gamma,
                           const std::vector<float>& mass_density)
{
    float m = 0.0f;
    Vector3f p_first = Vector3f::Zero();
    Matrix3f I_tensor = Matrix3f::Zero();

    for (int i = 0; i < (int)gamma.size(); i++) {
        float rho = mass_density[i];
        const Vector3f& g = gamma[i];

        m        += rho;
        p_first  += rho * g;
        I_tensor += rho * (g.squaredNorm() * Matrix3f::Identity() - g * g.transpose());
    }

    Matrix3f P = hat(p_first);

    Matrix6f M = Matrix6f::Zero();
    M.block<3,3>(0,0) = I_tensor;
    M.block<3,3>(0,3) = P;
    M.block<3,3>(3,0) = P.transpose();
    M.block<3,3>(3,3) = m * Matrix3f::Identity();
    return M;
}

// Algorithm 4: Fluid added mass matrix M_f
//
// Derived from calc_fluid_momentum: for rigid motion gamma_prime_ijk = omega x gamma_ijk + v,
// dot(gamma_prime, n) = dot(q, omega) + dot(n, v)  where q = gamma x n.
// Substituting gives M_f * [omega; v] = [l_f; p_f]:
//
//   top-left  (3x3): sum_f scale * q * q^T   (dl_f/domega)
//   top-right (3x3): sum_f scale * q * n^T   (dl_f/dv)
//   bot-left  (3x3): sum_f scale * n * q^T   (dp_f/domega)
//   bot-right (3x3): sum_f scale * n * n^T   (dp_f/dv)
//
// where scale = rho_f * delta * A_f.  M_f is symmetric.
Matrix6f calc_added_mass(const std::vector<Vector3f>& gamma,
                         const std::vector<Vector3i>& faces,
                         const std::vector<float>& face_areas,
                         const std::vector<Vector3f>& face_normals,
                         float fluid_density,
                         float delta)
{
    Matrix6f M = Matrix6f::Zero();

    for (int f = 0; f < (int)faces.size(); f++) {
        const int i = faces[f].x();
        const int j = faces[f].y();
        const int k = faces[f].z();

        const Vector3f centroid = (gamma[i] + gamma[j] + gamma[k]) / 3.0f;
        const Vector3f n = face_normals[f];
        const float    A = face_areas[f];
        const Vector3f q = centroid.cross(n);

        const float scale = fluid_density * delta * A;

        M.block<3,3>(0,0) += scale * q * q.transpose();
        M.block<3,3>(0,3) += scale * q * n.transpose();
        M.block<3,3>(3,0) += scale * n * q.transpose();
        M.block<3,3>(3,3) += scale * n * n.transpose();
    }

    return M;
}
