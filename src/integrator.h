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

    void LoadPose(const std::string& filepath, std::vector<Vector3f>& vertices, std::vector<Vector3i>& faces, std::vector<Edge>& edges);

    void CalculateFaceAttributes(std::vector<Vector3i>& faces, std::vector<Vector3f>& gamma);
    float CalculateDelta(std::vector<Vector3i>& faces, std::vector<Edge>& edges, std::vector<Vector3f>& gamma);

private:
    //std::vector<Vector3f> vertices;
    //std::vector<Vector3f> faces;

    std::vector<float> face_areas;
    std::vector<Vector3f> face_normals;
};

#endif // INTEGRATOR_H
