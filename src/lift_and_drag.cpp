#include <vector>
#include <Eigen/Dense>

using namespace Eigen;
using Vector6f = Eigen::Matrix<float, 6, 1>;

// Put two Vector3f halves into a Vector6f.
static inline Vector6f make_vec6(const Vector3f& rot, const Vector3f& trans) {
    Vector6f v;
    v << rot, trans;
    return v;
}

// =============================================================================
// External Forces
//
// Computes all external-force wrenches acting on the body as 6-vectors
// (torque, force) in the body frame.
//
//   // --- Mesh (body frame) ---
//   std::vector<Vector3f>  positions_k;  // vertex positions at time k   (body)
//   std::vector<Vector3f>  positions_k1; // vertex positions at time k+1 (body)
//   std::vector<Vector3i>  faces;        // triangle indices
//   std::vector<float>     mass_density; // per-vertex body mass density (kg/m^3)
//
//   // --- State ---
//   Matrix4f pose;                       // current pose (SE(3)): rotation A, translation b
//   Vector6f body_twist;                 // current body-frame twist (omega, v)
//
//   // --- Physical parameters ---
//   float fluid_density;                 // ambient fluid density (kg/m^3)
//   float timestep;                      // integration time step (s)
//   Vector3f gravity;                    // gravity vector in world space, e.g. (0,-9.81,0)
//   float total_mass;                    // effective mass: (rho_b - rho_f) * V
//
//   // --- Background flow ---
//   Vector3f background_flow;            // typically zero current
//
// =============================================================================


// Volume of a closed triangle mesh (body coordinates)
//
// Uses the divergence theorem:
//     V = (1/3) * sum_faces <gamma_ijk, n_ijk> * A_ijk
// where gamma_ijk is the face centroid and n_ijk is the outward unit normal.
float compute_volume(const std::vector<Vector3f>& positions,
                     const std::vector<Vector3i>& faces)
{
    float V = 0.0f;
    for (int f = 0; f < (int)faces.size(); f++) {
        const int i = faces[f](0);
        const int j = faces[f](1);
        const int k = faces[f](2);

        const Vector3f p_i = positions[i];
        const Vector3f p_j = positions[j];
        const Vector3f p_k = positions[k];

        const Vector3f c = (p_i + p_j + p_k) / 3.0f;       // centroid
        const Vector3f N = (p_j - p_i).cross(p_k - p_i);   // 2*A*n
        const float Af   = 0.5f * N.norm();

        if (Af < 1e-12f) continue;

        const Vector3f n = N / (2.0f * Af);                // unit normal
        V += c.dot(n) * Af / 3.0f;
    }
    return V > 0.0f ? V : 0.0f;
}


// Lift and Drag  (Alg. 2 / Prop. 1)
Vector6f calc_lift_and_drag(const std::vector<Vector3f>& positions_k,
                            const std::vector<Vector3f>& positions_k1,
                            const std::vector<Vector3i>& faces,
                            const Matrix4f& pose,
                            const Vector6f& body_twist,
                            float fluid_density,
                            float timestep,
                            const Vector3f& background_flow)
{
    // Pull the constant background flow into the body frame
    const Matrix3f A_T = pose.block<3,3>(0,0).transpose();
    const Vector3f u_bg_body = A_T * background_flow;

    // Unpack the body-frame twist body_twist = (omega, v)
    const Vector3f omega = body_twist.head<3>();
    const Vector3f v     = body_twist.tail<3>();

    // Accumulators for the total wrench in the body frame
    Vector3f torque = Vector3f::Zero();
    Vector3f force  = Vector3f::Zero();

    // Loop over every triangle of the surface mesh
    for (int f = 0; f < (int)faces.size(); f++) {
        const int i = faces[f](0);
        const int j = faces[f](1);
        const int k = faces[f](2);

        // Centroid of the triangle
        const Vector3f centroid = (positions_k[i] + positions_k[j] + positions_k[k]) / 3.0f;

        // Unnormalized normal = (p_j - p_i) x (p_k - p_i); area = 0.5 * |N|
        const Vector3f N_raw = (positions_k[j] - positions_k[i]).cross(positions_k[k] - positions_k[i]);
        const float A_ijk    = 0.5f * N_raw.norm();
        if (A_ijk < 1e-12f) continue;              
        const Vector3f n_ijk = N_raw / (2.0f * A_ijk);  

        // Face velocity in the body frame
        // (1) rigid contribution from the current twist Y
        const Vector3f v_rigid = omega.cross(centroid) + v;

        // (2) shape-change contribution from the pose sequence
        const Vector3f centroid_prime =
            ((positions_k1[i] - positions_k[i]) +
             (positions_k1[j] - positions_k[j]) +
             (positions_k1[k] - positions_k[k])) / (3.0f * timestep);

        const Vector3f v_face_body = v_rigid + centroid_prime;

        // Relative velocity of the surface w.r.t. the fluid
        const Vector3f u = v_face_body - u_bg_body;
        const float u_mag = u.norm();
        if (u_mag < 1e-12f) continue;             

        // Prop. 1:  F_face = -rho_f * A * |u| * <u,n> * n
        const float un = u.dot(n_ijk);
        const Vector3f f_face = -fluid_density * A_ijk * u_mag * un * n_ijk;

        // Accumulate as a body-frame wrench
        // Force applied at the centroid contributes torque centroid x f_face.
        force  += f_face;
        torque += centroid.cross(f_face);
    }

    return make_vec6(torque, force);
}


// A uniform downward Gravity force in world space, applied at the center of
// mass. Assumes the mesh's origin is already at the center of mass.
Vector6f calc_gravity(const Matrix4f& pose,
                      const Vector3f& gravity,
                      float total_mass)
{
    const Matrix3f A_T = pose.block<3,3>(0,0).transpose();

    const Vector3f f_world  = total_mass * gravity;
    const Vector3f f_body   = A_T * f_world;
    const Vector3f tau_body = Vector3f::Zero();

    return make_vec6(tau_body, f_body);
}


// Total external force, used from the Newton solve's residual evaluation.
Vector6f calc_total_force(const std::vector<Vector3f>& positions_k,
                          const std::vector<Vector3f>& positions_k1,
                          const std::vector<Vector3i>& faces,
                          const Matrix4f& pose,
                          const Vector6f& body_twist,
                          float fluid_density,
                          float timestep,
                          const Vector3f& gravity,
                          float total_mass,
                          const Vector3f& background_flow)
{
    // Initialize zero force
    Vector6f externalForce = Vector6f::Zero();

    // Calculate lift and drag from movement in the fluid
    externalForce += calc_lift_and_drag(positions_k, positions_k1, faces, pose, body_twist,
                                        fluid_density, timestep, background_flow);

    // Apply gravity (total_mass accounts for buoyancy)
    externalForce += calc_gravity(pose, gravity, total_mass);

    return externalForce;
}
