#ifndef INTEGRATOR_H
#define INTEGRATOR_H

using namespace Eigen;
using Vector6f = Eigen::Matrix<float, 6, 1>;

// Put two Vector3f halves into a Vector6f.
static inline Vector6f make_vec6(const Vector3f& rot, const Vector3f& trans) {
    Vector6f v;
    v << rot, trans;
    return v;
}

class Integrator
{
public:
    Integrator();
};

#endif // INTEGRATOR_H
