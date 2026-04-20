#ifndef INERTIA_H
#define INERTIA_H

#include "integrator.h"

// Algorithm 3: Time-varying body inertia matrix
// Returns the 6x6 symmetric matrix M_b such that M_b * [omega; v] = [l_b; p_b]
Matrix6f calc_body_inertia(const std::vector<Vector3f>& gamma,
                           const std::vector<float>& mass_density);

// Algorithm 4: Fluid added mass matrix
// Returns the 6x6 symmetric matrix M_f such that M_f * [omega; v] = [l_f; p_f]
Matrix6f calc_added_mass(const std::vector<Vector3f>& gamma,
                         const std::vector<Vector3i>& faces,
                         const std::vector<float>& face_areas,
                         const std::vector<Vector3f>& face_normals,
                         float fluid_density,
                         float delta);

#endif // INERTIA_H
