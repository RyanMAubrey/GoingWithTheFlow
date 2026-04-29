#ifndef SE3_H
#define SE3_H

#include "integrator.h"

// =============================================================================
// SE(3) / se(3) Utilities
//
// These implement the Cayley map and its differential as described in
// Appendix B of "Going with the Flow".
//
// Conventions:
//   * Vector6f layout: indices 0..2 = rotational (ω), 3..5 = translational (v)
//   * Matrix4f SE(3) element: top-left 3×3 = rotation A, col 3 top 3 = translation b
//   * hat(ω) is the 3×3 skew-symmetric matrix such that hat(ω) * x = ω × x
// =============================================================================

// Skew-symmetric (hat) matrix: hat(w) * x == w.cross(x)
inline Matrix3f hat(const Vector3f& w) {
    Matrix3f m;
    m <<     0, -w.z(),  w.y(),
         w.z(),      0, -w.x(),
        -w.y(),  w.x(),      0;
    return m;
}

// Cayley map  τ : se(3) → SE(3)
inline Matrix4f cayley_map(const Vector6f& Y) {
    const Vector3f w = Y.head<3>();
    const Vector3f v = Y.tail<3>();

    const Matrix3f W = hat(w);
    const float w_sq = w.squaredNorm();
    const float s = 4.0f / (4.0f + w_sq);

    // Rotation: A = I + s * (W + ½ W²)
    const Matrix3f A = Matrix3f::Identity() + s * (W + 0.5f * W * W);

    // Translation: b = (s/2) * (2v + ω×v)  =  s * v + (s/2) * (ω×v)
    const Vector3f b = s * v + (s * 0.5f) * w.cross(v);

    Matrix4f g = Matrix4f::Identity();
    g.block<3,3>(0,0) = A;
    g.block<3,1>(0,3) = b;
    return g;
}

// Inverse right-trivialized differential of the Cayley map, transposed.
inline Matrix6f dtau_inv_star(const Vector6f& Y) {
    const Vector3f w = Y.head<3>();
    const Vector3f v = Y.tail<3>();

    const Matrix3f W = hat(w);
    const Matrix3f V = hat(v);
    const Matrix3f I = Matrix3f::Identity();

    // Non-transposed blocks (Eqn. 17)
    const Matrix3f TL = I - 0.5f * W + 0.25f * (w * w.transpose());
    const Matrix3f TR = Matrix3f::Zero();
    const Matrix3f BL = -0.5f * (V - 0.5f * W * V);
    const Matrix3f BR = I - 0.5f * W;

    // Assemble the non-transposed 6×6
    Matrix6f D = Matrix6f::Zero();
    D.block<3,3>(0,0) = TL;
    D.block<3,3>(0,3) = TR;
    D.block<3,3>(3,0) = BL;
    D.block<3,3>(3,3) = BR;

    // Return the transpose
    return D.transpose();
}

#endif // SE3_H
