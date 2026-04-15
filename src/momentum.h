#ifndef MOMENTUM_H
#define MOMENTUM_H

#include <vector>
#include <Eigen/Dense>
#include <integrator.h>

class Momentum
{
public:
    Momentum();

    Vector6f calc_body_momentum(std::vector<Vector3f>& gamma_k, std::vector<Vector3f>& gamma_k1,
                                std::vector<float>& mass_density, float h);
    Vector6f calc_fluid_momentum(std::vector<Vector3f>& gamma_k, std::vector<Vector3f>& gamma_k1, std::vector<Vector3i>& faces,
                                 std::vector<float> face_areas, std::vector<Vector3f> face_normals,
                                 std::vector<Edge>& edges, float rho_f, float h, float delta);
};

#endif // MOMENTUM_H
