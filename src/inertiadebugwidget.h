#pragma once

#include <QOpenGLWidget>
#include <QOpenGLFunctions>
#include <QMatrix4x4>
#include <Eigen/Dense>

#include "mesh_loader.h"
#include "inertia.h"
#include "integrator.h"

// Visualizes the principal axes of the 3x3 sub-blocks of the body inertia
// tensor M_b and added mass tensor M_f.
//
// Eigenvectors of a symmetric 3x3 inertia block are the body's "natural"
// rotation/translation axes; eigenvalues are how much momentum each axis
// produces under unit motion. Drawn as arrows from the body center, with
// length proportional to the eigenvalue (normalized per tensor).
//
//   green  = M_b top-left   (body rotational inertia)
//   yellow = M_f top-left   (added rotational inertia)
//   cyan   = M_f bottom-right (added translational mass)
//
// (M_b bottom-right is m*I3 -- isotropic, no axes worth drawing.)
class InertiaDebugWidget : public QOpenGLWidget, protected QOpenGLFunctions
{
    Q_OBJECT

public:
    explicit InertiaDebugWidget(int frame = 1, QWidget *parent = nullptr);

protected:
    void initializeGL() override;
    void resizeGL(int w, int h) override;
    void paintGL() override;

private:
    int m_frame = 1;
    TriMesh m_mesh;
    QMatrix4x4 m_proj;

    Eigen::Vector3f m_center = Eigen::Vector3f::Zero();
    float m_scale = 1.0f;

    Eigen::Matrix3f m_body_rot_axes   = Eigen::Matrix3f::Identity();
    Eigen::Vector3f m_body_rot_vals   = Eigen::Vector3f::Zero();
    Eigen::Matrix3f m_fluid_rot_axes  = Eigen::Matrix3f::Identity();
    Eigen::Vector3f m_fluid_rot_vals  = Eigen::Vector3f::Zero();
    Eigen::Matrix3f m_fluid_trans_axes = Eigen::Matrix3f::Identity();
    Eigen::Vector3f m_fluid_trans_vals = Eigen::Vector3f::Zero();

    void loadFrames();
    void computeViewTransform();
    Eigen::Vector3f toViewSpace(const Eigen::Vector3f& v) const;

    void drawMesh(const TriMesh& mesh);
    void drawPrincipalAxes(const Eigen::Matrix3f& axes,
                           const Eigen::Vector3f& vals,
                           float r, float g, float b);
};
