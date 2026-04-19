#ifndef INTEGRATOR_H
#define INTEGRATOR_H

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

struct Edge {
    int i;
    int j;
    float alpha; // [-pi, pi]
};

class Integrator
{
public:
    Integrator();

    void Simulate();

    void LoadPose(const std::string& filepath, std::vector<Vector3f>& vertices, std::vector<Vector3i>& faces, std::vector<Edge>& edges,
                  std::vector<float>& face_areas, std::vector<Vector3f>& face_normals);

    void CalculateFaceAttributes(std::vector<Vector3i>& faces, std::vector<Vector3f>& gamma,
                                 std::vector<float>& face_areas, std::vector<Vector3f>& face_normals);
    float CalculateDelta(std::vector<Vector3i>& faces, std::vector<Edge>& edges, std::vector<Vector3f>& gamma, std::vector<float>& face_areas);

private:
    int total_poses = 6;
    float h = 7.0; // 8 frames per pose (may need to change)
    float rho_f = 1.0f; // Fluid density (may change)

    bool HasEdge(Vector3i& face, int i, int j);
};

#endif // INTEGRATOR_H
