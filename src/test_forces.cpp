#include "integrator.h"
#include "inertia.h"
#include "mesh_loader.h"
#include <iostream>
#include <cassert>
#include <cmath>

// Declarations
float compute_volume(const std::vector<Vector3f>& pos_k,
                     const std::vector<Vector3i>& faces);

Vector6f calc_lift_and_drag(const std::vector<Vector3f>& pos_k,
                            const std::vector<Vector3f>& pos_k1,
                            const std::vector<Vector3i>& faces,
                            const Matrix4f& pose,
                            const Vector6f& body_twist,
                            float fluid_density,
                            float timestep,
                            const Vector3f& background_flow);

Vector6f calc_gravity(const Matrix4f& pose,
                      const Vector3f& gravity,
                      float total_mass);

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

// ---- Helpers ----------------------------------------------------------------

static bool approx(float a, float b, float tol = 1e-6f) {
    return std::abs(a - b) < tol;
}

static bool approx_vec(const Vector3f& a, const Vector3f& b, float tol = 1e-6f) {
    return (a - b).norm() < tol;
}

static bool approx_vec6(const Vector6f& a, const Vector6f& b, float tol = 1e-6f) {
    return (a - b).norm() < tol;
}

static void print_vec6(const std::string& label, const Vector6f& v) {
    std::cout << "  " << label << ": ("
              << v[0] << ", " << v[1] << ", " << v[2] << " | "
              << v[3] << ", " << v[4] << ", " << v[5] << ")" << std::endl;
}

// A 1x1 flat plate centered at the origin in the XY plane, normal = +Z.
static void make_unit_plate(std::vector<Vector3f>& verts, std::vector<Vector3i>& faces) {
    verts = {
        Vector3f(-0.5f, -0.5f, 0.0f),
        Vector3f( 0.5f, -0.5f, 0.0f),
        Vector3f( 0.5f,  0.5f, 0.0f),
        Vector3f(-0.5f,  0.5f, 0.0f)
    };
    faces = {
        Vector3i(0, 1, 2),
        Vector3i(0, 2, 3)
    };
}

// A unit cube centered at the origin, closed mesh with outward normals.
// 8 vertices, 12 triangles.
static void make_unit_cube(std::vector<Vector3f>& verts, std::vector<Vector3i>& faces) {
    float s = 0.5f;
    verts = {
        Vector3f(-s, -s, -s),  // 0
        Vector3f( s, -s, -s),  // 1
        Vector3f( s,  s, -s),  // 2
        Vector3f(-s,  s, -s),  // 3
        Vector3f(-s, -s,  s),  // 4
        Vector3f( s, -s,  s),  // 5
        Vector3f( s,  s,  s),  // 6
        Vector3f(-s,  s,  s)   // 7
    };
    faces = {
        // -Z face
        Vector3i(0, 3, 2), Vector3i(0, 2, 1),
        // +Z face
        Vector3i(4, 5, 6), Vector3i(4, 6, 7),
        // -X face
        Vector3i(0, 4, 7), Vector3i(0, 7, 3),
        // +X face
        Vector3i(1, 2, 6), Vector3i(1, 6, 5),
        // -Y face
        Vector3i(0, 1, 5), Vector3i(0, 5, 4),
        // +Y face
        Vector3i(2, 3, 7), Vector3i(2, 7, 6)
    };
}


// ---- Test group 1: zero-input sanity ----------------------------------------

void test_zero_twist_zero_force() {
    std::cout << "test_zero_twist_zero_force... ";

    std::vector<Vector3f> verts;
    std::vector<Vector3i> faces;
    make_unit_plate(verts, faces);

    Matrix4f pose = Matrix4f::Identity();
    Vector6f twist = Vector6f::Zero();
    float rho_f = 998.0f;
    float h = 1.0f / 200.0f;
    Vector3f bg = Vector3f::Zero();

    // Rigid body: pos_k1 = pos_k, twist = 0, no background flow
    // => every face has u = 0 => force must be zero
    Vector6f F = calc_lift_and_drag(verts, verts, faces, pose, twist, rho_f, h, bg);

    assert(approx_vec6(F, Vector6f::Zero()));
    std::cout << "PASSED" << std::endl;
}

void test_gravity_identity_pose() {
    std::cout << "test_gravity_identity_pose... ";

    Matrix4f pose = Matrix4f::Identity();
    Vector3f g_world(0.0f, -9.81f, 0.0f);
    float mass = 2.0f;

    Vector6f F = calc_gravity(pose, g_world, mass);

    // At identity pose, body frame = world frame
    // f_body should be (0, -19.62, 0), torque = 0
    Vector3f expected_f(0.0f, -19.62f, 0.0f);
    assert(approx_vec(F.tail<3>(), expected_f, 1e-4f));
    assert(approx_vec(F.head<3>(), Vector3f::Zero()));
    std::cout << "PASSED" << std::endl;
}

void test_gravity_rotated_pose() {
    std::cout << "test_gravity_rotated_pose... ";

    // Rotate body 90 degrees about the X axis:
    // world Y -> body Z, world Z -> body -Y
    Matrix4f pose = Matrix4f::Identity();
    float angle = M_PI / 2.0f;
    // Rotation about X:
    //   [1   0     0  ]
    //   [0  cos  -sin ]
    //   [0  sin   cos ]
    pose(1,1) =  std::cos(angle);
    pose(1,2) = -std::sin(angle);
    pose(2,1) =  std::sin(angle);
    pose(2,2) =  std::cos(angle);

    Vector3f g_world(0.0f, -9.81f, 0.0f);
    float mass = 1.0f;

    Vector6f F = calc_gravity(pose, g_world, mass);

    // gravity = (0, -9.81, 0) in world
    // A^T * gravity: with 90-deg X rotation, A^T maps world (0,1,0) to body (0,0,1)
    // so world (0,-9.81,0) maps to body (0, 0, +9.81)
    Vector3f expected_f(0.0f, 0.0f, 9.81f);
    assert(approx_vec(F.tail<3>(), expected_f, 1e-4f));
    assert(approx_vec(F.head<3>(), Vector3f::Zero()));

    print_vec6("gravity (rotated)", F);
    std::cout << "PASSED" << std::endl;
}

void test_gravity_zero_mass() {
    std::cout << "test_gravity_zero_mass... ";

    // Neutrally buoyant: total_mass = (rho_b - rho_f) * V = 0
    Matrix4f pose = Matrix4f::Identity();
    Vector3f g_world(0.0f, -9.81f, 0.0f);

    Vector6f F = calc_gravity(pose, g_world, 0.0f);

    assert(approx_vec6(F, Vector6f::Zero()));
    std::cout << "PASSED" << std::endl;
}


// ---- Test group 2: flat plate analytical ------------------------------------

void test_plate_normal_motion() {
    std::cout << "test_plate_normal_motion... ";

    // 1x1 plate in XY plane, normal = +Z
    // Push it in -Z direction (into the fluid, normal-first)
    std::vector<Vector3f> verts;
    std::vector<Vector3i> faces;
    make_unit_plate(verts, faces);

    Matrix4f pose = Matrix4f::Identity();
    Vector6f twist = Vector6f::Zero();
    twist[5] = -1.0f;  // v_z = -1 m/s (moving in -Z)

    float rho_f = 1.0f;  // unit density for easy math
    float h = 1.0f / 200.0f;
    Vector3f bg = Vector3f::Zero();

    Vector6f F = calc_lift_and_drag(verts, verts, faces, pose, twist, rho_f, h, bg);

    // By Prop. 1: F_face = -rho_f * A * |u| * <u,n> * n
    //
    // u = (0, 0, -1), n = (0, 0, 1) for both triangles
    // |u| = 1, <u,n> = -1, A_total = 1.0
    //
    // F = -1.0 * 1.0 * 1.0 * (-1.0) * (0,0,1) = (0, 0, 1)
    // (force opposes the motion — pushes back in +Z)
    //
    // Torque: centroid = (0,0,0) for the whole plate => torque = 0

    Vector3f expected_f(0.0f, 0.0f, 1.0f);
    assert(approx_vec(F.tail<3>(), expected_f, 1e-4f));
    assert(approx_vec(F.head<3>(), Vector3f::Zero(), 1e-4f));

    print_vec6("plate normal motion", F);
    std::cout << "PASSED" << std::endl;
}

void test_plate_tangential_motion() {
    std::cout << "test_plate_tangential_motion... ";

    // Same plate, but slide it in +X (tangent to the surface)
    std::vector<Vector3f> verts;
    std::vector<Vector3i> faces;
    make_unit_plate(verts, faces);

    Matrix4f pose = Matrix4f::Identity();
    Vector6f twist = Vector6f::Zero();
    twist[3] = 1.0f;  // v_x = +1 m/s

    float rho_f = 998.0f;
    float h = 1.0f / 200.0f;
    Vector3f bg = Vector3f::Zero();

    Vector6f F = calc_lift_and_drag(verts, verts, faces, pose, twist, rho_f, h, bg);

    // u = (1,0,0), n = (0,0,1)
    // <u,n> = 0 => F_face = 0 for every triangle
    // Tangential motion produces ZERO lift+drag (no viscous skin friction in this model)

    assert(approx_vec6(F, Vector6f::Zero(), 1e-4f));
    std::cout << "PASSED" << std::endl;
}

void test_plate_45_degree_motion() {
    std::cout << "test_plate_45_degree_motion... ";

    // Plate in XY plane, moving at 45 degrees between X and Z
    std::vector<Vector3f> verts;
    std::vector<Vector3i> faces;
    make_unit_plate(verts, faces);

    Matrix4f pose = Matrix4f::Identity();
    Vector6f twist = Vector6f::Zero();
    float speed = 1.0f;
    twist[3] = speed * std::cos(M_PI / 4.0f);  // v_x
    twist[5] = -speed * std::sin(M_PI / 4.0f); // v_z (into -Z)

    float rho_f = 1.0f;
    float h = 1.0f / 200.0f;
    Vector3f bg = Vector3f::Zero();

    Vector6f F = calc_lift_and_drag(verts, verts, faces, pose, twist, rho_f, h, bg);

    // u = (cos45, 0, -sin45), n = (0, 0, 1)
    // |u| = 1.0, <u,n> = -sin45
    // F = -1.0 * 1.0 * 1.0 * (-sin45) * (0,0,1) = (0, 0, sin45)
    //
    // Key check: force is purely in the normal direction (Z), not along u.
    // This is the paper's Prop. 1 in action — the combined lift+drag resolves
    // to a single force along the face normal.

    float sin45 = std::sin(M_PI / 4.0f);
    Vector3f expected_f(0.0f, 0.0f, sin45);

    assert(approx_vec(F.tail<3>(), expected_f, 1e-4f));
    assert(approx_vec(F.head<3>(), Vector3f::Zero(), 1e-4f));

    print_vec6("plate 45-deg motion", F);
    std::cout << "PASSED" << std::endl;
}

void test_plate_background_flow() {
    std::cout << "test_plate_background_flow... ";

    // Stationary plate, but background flow pushes fluid in -Z
    // Should be equivalent to moving the plate in +Z through still fluid
    std::vector<Vector3f> verts;
    std::vector<Vector3i> faces;
    make_unit_plate(verts, faces);

    Matrix4f pose = Matrix4f::Identity();
    Vector6f twist = Vector6f::Zero();  // plate is stationary

    float rho_f = 1.0f;
    float h = 1.0f / 200.0f;
    Vector3f bg(0.0f, 0.0f, -1.0f);  // flow in -Z

    Vector6f F = calc_lift_and_drag(verts, verts, faces, pose, twist, rho_f, h, bg);

    // u = v_face - u_bg = 0 - (0,0,-1) = (0, 0, 1)
    // n = (0, 0, 1), <u,n> = 1, |u| = 1
    // F = -1.0 * 1.0 * 1.0 * 1.0 * (0,0,1) = (0, 0, -1)
    //
    // Wind blowing in -Z pushes plate in -Z.

    Vector3f expected_f(0.0f, 0.0f, -1.0f);
    assert(approx_vec(F.tail<3>(), expected_f, 1e-4f));

    print_vec6("plate bg flow", F);
    std::cout << "PASSED" << std::endl;
}


// ---- Test group 3: closed mesh / volume / symmetry --------------------------

void test_unit_cube_volume() {
    std::cout << "test_unit_cube_volume... ";

    std::vector<Vector3f> verts;
    std::vector<Vector3i> faces;
    make_unit_cube(verts, faces);

    float V = compute_volume(verts, faces);

    // Unit cube: volume should be 1.0
    assert(approx(V, 1.0f, 1e-4f));
    std::cout << "PASSED (V = " << V << ")" << std::endl;
}

void test_cube_stationary_symmetric() {
    std::cout << "test_cube_stationary_symmetric... ";

    // A stationary cube in a uniform flow should get a net force
    // along the flow direction (by symmetry, no transverse force)
    std::vector<Vector3f> verts;
    std::vector<Vector3i> faces;
    make_unit_cube(verts, faces);

    Matrix4f pose = Matrix4f::Identity();
    Vector6f twist = Vector6f::Zero();
    float rho_f = 1.0f;
    float h = 1.0f / 200.0f;
    Vector3f bg(1.0f, 0.0f, 0.0f);  // flow in +X

    Vector6f F = calc_lift_and_drag(verts, verts, faces, pose, twist, rho_f, h, bg);

    // By symmetry about the XZ and XY planes:
    //   f_y should be ~0
    //   f_z should be ~0
    //   f_x should be negative (drag opposes flow, pushes cube in -X? no —
    //         the cube is stationary and the flow is moving in +X, so
    //         u = -bg = (-1,0,0). On the -X face, <u,n> = <(-1,0,0),(-1,0,0)> = 1,
    //         F = -rho*A*|u|*1*(-1,0,0) = (1,0,0). On the +X face,
    //         <u,n> = <(-1,0,0),(1,0,0)> = -1, F = -rho*A*1*(-1)*(1,0,0) = (1,0,0).
    //         Both faces push in +X. This is correct: a headwind pushes the object
    //         downwind.
    //   torque should be ~0 (symmetric cube centered at origin)

    assert(std::abs(F.tail<3>().y()) < 1e-4f);  // f_y ~ 0
    assert(std::abs(F.tail<3>().z()) < 1e-4f);  // f_z ~ 0
    assert(std::abs(F.head<3>().norm()) < 1e-4f); // torque ~ 0 (symmetric)
    assert(F.tail<3>().x() > 0.0f);  // net force in +X (object pushed downwind)

    print_vec6("cube in +X flow", F);
    std::cout << "PASSED" << std::endl;
}

void test_shape_change_produces_thrust() {
    std::cout << "test_shape_change_produces_thrust... ";

    // A plate whose vertices move in +Z between k and k+1
    // should experience a force in -Z (reaction to pushing fluid)
    std::vector<Vector3f> verts_k, verts_k1;
    std::vector<Vector3i> faces;
    make_unit_plate(verts_k, faces);

    // Move all vertices 0.1 in +Z for the next frame
    verts_k1 = verts_k;
    for (auto& v : verts_k1) {
        v.z() += 0.1f;
    }

    Matrix4f pose = Matrix4f::Identity();
    Vector6f twist = Vector6f::Zero();  // no rigid motion yet
    float rho_f = 1.0f;
    float h = 1.0f / 200.0f;
    Vector3f bg = Vector3f::Zero();

    Vector6f F = calc_lift_and_drag(verts_k, verts_k1, faces, pose, twist, rho_f, h, bg);

    // centroid_prime = (0, 0, 0.1) / (h) ... well, per vertex the diff is 0.1,
    // so centroid_prime = (0, 0, 0.1/h) ... wait: the formula divides by (3*h)
    // but sums 3 vertex diffs, so it's (3 * 0.1) / (3 * h) = 0.1 / h.
    //
    // u = v_rigid + centroid_prime = (0,0,0) + (0, 0, 0.1/h)
    // n = (0, 0, 1), <u,n> = 0.1/h, |u| = 0.1/h
    //
    // F = -rho * A * |u| * <u,n> * n = -1 * 1 * (0.1/h) * (0.1/h) * (0,0,1)
    //   = (0, 0, -(0.1/h)^2)
    //
    // Key check: force is in -Z (opposing the shape change direction)

    assert(F.tail<3>().z() < 0.0f);

    float expected_fz = -(0.1f / h) * (0.1f / h);
    assert(approx(F.tail<3>().z(), expected_fz, std::abs(expected_fz) * 1e-3f));

    print_vec6("shape-change thrust", F);
    std::cout << "PASSED" << std::endl;
}


// ---- Test group 4: total force integration ----------------------------------

void test_total_force_neutrally_buoyant() {
    std::cout << "test_total_force_neutrally_buoyant... ";

    // Neutrally buoyant stationary rigid body in quiescent fluid
    // => every force is zero
    std::vector<Vector3f> verts;
    std::vector<Vector3i> faces;
    make_unit_cube(verts, faces);

    Matrix4f pose = Matrix4f::Identity();
    Vector6f twist = Vector6f::Zero();
    float rho_f = 998.0f;
    float h = 1.0f / 200.0f;
    Vector3f gravity(0.0f, -9.81f, 0.0f);
    float total_mass = 0.0f;  // neutrally buoyant
    Vector3f bg = Vector3f::Zero();

    Vector6f F = calc_total_force(verts, verts, faces, pose, twist,
                                   rho_f, h, gravity, total_mass, bg);

    assert(approx_vec6(F, Vector6f::Zero()));
    std::cout << "PASSED" << std::endl;
}


// ---- Test group 5: inertia matrices (Alg. 3 & 4) ---------------------------

void test_body_inertia_symmetry() {
    std::cout << "test_body_inertia_symmetry... ";

    std::vector<Vector3f> verts;
    std::vector<Vector3i> faces;
    make_unit_cube(verts, faces);

    std::vector<float> rho(verts.size(), 1.0f);
    Matrix6f M = calc_body_inertia(verts, rho);

    assert((M - M.transpose()).norm() < 1e-5f);
    std::cout << "PASSED" << std::endl;
}

void test_body_inertia_matches_momentum() {
    std::cout << "test_body_inertia_matches_momentum... ";

    // For pure rigid motion, M_b * twist should equal calc_body_momentum result.
    std::vector<Vector3f> verts_k;
    std::vector<Vector3i> faces;
    make_unit_plate(verts_k, faces);

    // Rigid twist: spin about Z and translate in X
    Vector6f twist = Vector6f::Zero();
    twist[2] = 0.5f;  // omega_z
    twist[3] = 1.0f;  // v_x

    float h = 1.0f / 200.0f;
    std::vector<float> rho(verts_k.size(), 2.0f);

    // Build k+1 by applying rigid twist for one step
    std::vector<Vector3f> verts_k1 = verts_k;
    Vector3f omega = twist.head<3>();
    Vector3f v     = twist.tail<3>();
    for (int i = 0; i < (int)verts_k.size(); i++) {
        verts_k1[i] = verts_k[i] + h * (omega.cross(verts_k[i]) + v);
    }

    // Direct momentum computation (Briana's formula)
    Vector3f l_direct = Vector3f::Zero(), p_direct = Vector3f::Zero();
    for (int i = 0; i < (int)verts_k.size(); i++) {
        Vector3f gp = (verts_k1[i] - verts_k[i]) / h;
        l_direct += rho[i] * verts_k[i].cross(gp);
        p_direct += rho[i] * gp;
    }
    Vector6f mom_direct = make_vec6(l_direct, p_direct);

    // Matrix multiply
    Matrix6f M = calc_body_inertia(verts_k, rho);
    Vector6f mom_matrix = M * twist;

    assert(approx_vec6(mom_direct, mom_matrix, 1e-4f));
    std::cout << "PASSED" << std::endl;
}

void test_added_mass_symmetry() {
    std::cout << "test_added_mass_symmetry... ";

    std::vector<Vector3f> verts;
    std::vector<Vector3i> faces;
    make_unit_cube(verts, faces);

    // Compute face attributes inline (mirrors Integrator::CalculateFaceAttributes)
    std::vector<float> areas(faces.size());
    std::vector<Vector3f> normals(faces.size());
    for (int f = 0; f < (int)faces.size(); f++) {
        Vector3f e1 = verts[faces[f].y()] - verts[faces[f].x()];
        Vector3f e2 = verts[faces[f].z()] - verts[faces[f].x()];
        Vector3f N  = e1.cross(e2);
        areas[f]   = 0.5f * N.norm();
        normals[f] = N.normalized();
    }

    float delta = 1.0f;
    Matrix6f M = calc_added_mass(verts, faces, areas, normals, 1.0f, delta);

    assert((M - M.transpose()).norm() < 1e-5f);
    std::cout << "PASSED" << std::endl;
}

void test_added_mass_matches_fluid_momentum() {
    std::cout << "test_added_mass_matches_fluid_momentum... ";

    // For pure rigid motion, M_f * twist should equal calc_fluid_momentum result.
    std::vector<Vector3f> verts_k;
    std::vector<Vector3i> faces;
    make_unit_plate(verts_k, faces);

    std::vector<float> areas(faces.size());
    std::vector<Vector3f> normals(faces.size());
    for (int f = 0; f < (int)faces.size(); f++) {
        Vector3f e1 = verts_k[faces[f].y()] - verts_k[faces[f].x()];
        Vector3f e2 = verts_k[faces[f].z()] - verts_k[faces[f].x()];
        Vector3f N  = e1.cross(e2);
        areas[f]   = 0.5f * N.norm();
        normals[f] = N.normalized();
    }

    // Rigid twist: translate in Z (normal to plate)
    Vector6f twist = Vector6f::Zero();
    twist[5] = 1.0f;  // v_z

    float h     = 1.0f / 200.0f;
    float rho_f = 1.0f;
    float delta = 1.0f;

    // Build k+1 by applying rigid twist
    std::vector<Vector3f> verts_k1 = verts_k;
    Vector3f omega = twist.head<3>();
    Vector3f v_t   = twist.tail<3>();
    for (int i = 0; i < (int)verts_k.size(); i++) {
        verts_k1[i] = verts_k[i] + h * (omega.cross(verts_k[i]) + v_t);
    }

    // Direct fluid momentum
    Vector3f l_direct = Vector3f::Zero(), p_direct = Vector3f::Zero();
    for (int f = 0; f < (int)faces.size(); f++) {
        int i = faces[f].x(), j = faces[f].y(), k = faces[f].z();
        Vector3f g  = (verts_k[i]  + verts_k[j]  + verts_k[k])  / 3.0f;
        Vector3f gp = ((verts_k1[i]-verts_k[i]) + (verts_k1[j]-verts_k[j]) + (verts_k1[k]-verts_k[k])) / (3.0f * h);
        float un = gp.dot(normals[f]);
        l_direct += rho_f * delta * un * g.cross(normals[f]) * areas[f];
        p_direct += rho_f * delta * un * normals[f] * areas[f];
    }
    Vector6f mom_direct = make_vec6(l_direct, p_direct);

    // Matrix multiply
    Matrix6f M = calc_added_mass(verts_k, faces, areas, normals, rho_f, delta);
    Vector6f mom_matrix = M * twist;

    assert(approx_vec6(mom_direct, mom_matrix, 1e-4f));
    std::cout << "PASSED" << std::endl;
}


// ---- main -------------------------------------------------------------------

int main() {
    std::cout << "=== Forces Unit Tests ===" << std::endl;
    std::cout << std::endl;

    std::cout << "--- Group 1: Zero-Input ---" << std::endl;
    test_zero_twist_zero_force();
    test_gravity_identity_pose();
    test_gravity_rotated_pose();
    test_gravity_zero_mass();
    std::cout << std::endl;

    std::cout << "--- Group 2: Flat Plate ---" << std::endl;
    test_plate_normal_motion();
    test_plate_tangential_motion();
    test_plate_45_degree_motion();
    test_plate_background_flow();
    std::cout << std::endl;

    std::cout << "--- Group 3: Closed Mesh / Volume / Symmetry ---" << std::endl;
    test_unit_cube_volume();
    test_cube_stationary_symmetric();
    test_shape_change_produces_thrust();
    std::cout << std::endl;

    std::cout << "--- Group 4: Total Force Integration ---" << std::endl;
    test_total_force_neutrally_buoyant();
    std::cout << std::endl;

    std::cout << "--- Group 5: Inertia Matrices (Alg. 3 & 4) ---" << std::endl;
    test_body_inertia_symmetry();
    test_body_inertia_matches_momentum();
    test_added_mass_symmetry();
    test_added_mass_matches_fluid_momentum();
    std::cout << std::endl;

    std::cout << "=== ALL TESTS PASSED ===" << std::endl;
    return 0;
}
