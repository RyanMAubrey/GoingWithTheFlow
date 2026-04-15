#ifndef MOMENTUM_H
#define MOMENTUM_H

#include <vector>
#include <Eigen/Dense>

using namespace Eigen;
using Vector6f = Eigen::Matrix<float, 6, 1>;

struct Edge {
    int i;
    int j;
    float alpha; // [-pi, pi]
};

class Momentum
{
public:
    Momentum();
    Vector6f calc_body_momentum(std::vector<Vector3f>& gamma_k, std::vector<Vector3f>& gamma_k1,
                                std::vector<float>& mass_density, float h);
    Vector6f calc_fluid_momentum(std::vector<Vector3f>& gamma_k, std::vector<Vector3f>& gamma_k1,
                                 std::vector<Vector3i>& faces, std::vector<Edge>& edges,
                                 float rho_f, float h);
};

#endif // MOMENTUM_H
