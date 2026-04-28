#include "inertiadebugwidget.h"

#include <QDebug>
#include <QDir>
#include <QCoreApplication>
#include <algorithm>
#include <limits>
#include <vector>

#ifdef __APPLE__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif

using namespace Eigen;

InertiaDebugWidget::InertiaDebugWidget(int frame, QWidget *parent)
    : QOpenGLWidget(parent), m_frame(frame)
{
}

void InertiaDebugWidget::loadFrames() {
    qDebug() << "Current working dir:" << QDir::currentPath();
    qDebug() << "App dir:" << QCoreApplication::applicationDirPath();

    const std::string path = "turtle_poses/frame_" + std::to_string(m_frame) + ".obj";
    m_mesh = load_obj(path);

    qDebug() << "frame" << m_frame
             << "verts:" << m_mesh.vertices.size()
             << "faces:" << m_mesh.faces.size();

    // LoadPose recomputes the same vertices but also gives us face attributes
    // and edges, which calc_added_mass and CalculateDelta need.
    Integrator integrator;
    std::vector<Vector3f> vertices;
    std::vector<Vector3i> faces;
    std::vector<Edge>     edges;
    std::vector<float>    face_areas;
    std::vector<Vector3f> face_normals;

    integrator.LoadPose(path,
                        vertices, faces, edges, face_areas, face_normals);

    // Body inertia M_b (Algo 3)
    std::vector<float> mass_density(vertices.size(), 1.0f);
    Matrix6f M_b = calc_body_inertia(vertices, mass_density);

    // Added mass M_f (Algo 4 / paper Algo 5)
    float rho_f = 1.0f;
    float delta = integrator.CalculateDelta(faces, edges, vertices, face_areas);
    Matrix6f M_f = calc_added_mass(vertices, faces, face_areas, face_normals,
                                   rho_f, delta);

    qDebug() << "delta:" << delta;

    auto eigen3 = [](const Matrix3f& M, Matrix3f& axes, Vector3f& vals) {
        SelfAdjointEigenSolver<Matrix3f> es(M);
        axes = es.eigenvectors();
        vals = es.eigenvalues();
    };

    eigen3(M_b.block<3,3>(0,0), m_body_rot_axes,    m_body_rot_vals);
    eigen3(M_f.block<3,3>(0,0), m_fluid_rot_axes,   m_fluid_rot_vals);
    eigen3(M_f.block<3,3>(3,3), m_fluid_trans_axes, m_fluid_trans_vals);

    qDebug() << "M_b rotational eigenvalues:"
             << m_body_rot_vals[0] << m_body_rot_vals[1] << m_body_rot_vals[2];
    qDebug() << "M_f rotational eigenvalues:"
             << m_fluid_rot_vals[0] << m_fluid_rot_vals[1] << m_fluid_rot_vals[2];
    qDebug() << "M_f translational eigenvalues:"
             << m_fluid_trans_vals[0] << m_fluid_trans_vals[1] << m_fluid_trans_vals[2];

    computeViewTransform();
}

void InertiaDebugWidget::computeViewTransform() {
    if (m_mesh.vertices.empty()) {
        m_center = Vector3f::Zero();
        m_scale = 1.0f;
        return;
    }

    Vector3f minV( std::numeric_limits<float>::max(),
                   std::numeric_limits<float>::max(),
                   std::numeric_limits<float>::max() );
    Vector3f maxV(-std::numeric_limits<float>::max(),
                  -std::numeric_limits<float>::max(),
                  -std::numeric_limits<float>::max());

    for (const auto& v : m_mesh.vertices) {
        minV = minV.cwiseMin(v);
        maxV = maxV.cwiseMax(v);
    }

    m_center = 0.5f * (minV + maxV);

    Vector3f extent = maxV - minV;
    float maxExtent = std::max(extent.x(), std::max(extent.y(), extent.z()));

    m_scale = (maxExtent > 0.0f) ? 1.8f / maxExtent : 1.0f;

    qDebug() << "center:" << m_center.x() << m_center.y() << m_center.z();
    qDebug() << "scale:" << m_scale;
}

Vector3f InertiaDebugWidget::toViewSpace(const Vector3f& v) const {
    return (v - m_center) * m_scale;
}

void InertiaDebugWidget::initializeGL() {
    initializeOpenGLFunctions();

    glEnable(GL_DEPTH_TEST);
    glDisable(GL_CULL_FACE);
    glClearColor(0.1f, 0.1f, 0.1f, 1.0f);

    loadFrames();
}

void InertiaDebugWidget::resizeGL(int w, int h) {
    glViewport(0, 0, w, h);

    m_proj.setToIdentity();
    m_proj.perspective(45.0f, float(w) / float(std::max(1, h)), 0.1f, 100.0f);
}

void InertiaDebugWidget::drawMesh(const TriMesh& mesh) {
    glColor3f(0.75f, 0.75f, 0.75f);

    glBegin(GL_TRIANGLES);
    for (const auto& f : mesh.faces) {
        const auto a = toViewSpace(mesh.vertices[f[0]]);
        const auto b = toViewSpace(mesh.vertices[f[1]]);
        const auto c = toViewSpace(mesh.vertices[f[2]]);
        glVertex3f(a.x(), a.y(), a.z());
        glVertex3f(b.x(), b.y(), b.z());
        glVertex3f(c.x(), c.y(), c.z());
    }
    glEnd();
}

void InertiaDebugWidget::drawPrincipalAxes(const Matrix3f& axes,
                                           const Vector3f& vals,
                                           float r, float g, float b) {
    Vector3f center = toViewSpace(m_center);

    float max_val = vals.cwiseAbs().maxCoeff();
    if (max_val < 1e-12f) return;

    glColor3f(r, g, b);
    glLineWidth(4.0f);

    glBegin(GL_LINES);
    for (int k = 0; k < 3; k++) {
        Vector3f dir = axes.col(k).normalized();

        // Eigenvalue determines visual length; principal axes are unsigned lines,
        // so draw both halves through the body center.
        float vis_len = (std::abs(vals[k]) / max_val) * 0.8f;
        if (vis_len < 1e-4f) continue;

        Vector3f end_pos = center + dir * vis_len;
        Vector3f end_neg = center - dir * vis_len;

        // Main shaft (positive half + negative half)
        glVertex3f(center.x(), center.y(), center.z());
        glVertex3f(end_pos.x(), end_pos.y(), end_pos.z());
        glVertex3f(center.x(), center.y(), center.z());
        glVertex3f(end_neg.x(), end_neg.y(), end_neg.z());

        // Single arrowhead on the positive side just so the eigenvector
        // orientation is unambiguous in the rendering.
        Vector3f perp = dir.unitOrthogonal();
        float head = std::clamp(0.15f * vis_len, 0.02f, 0.06f);

        Vector3f left  = end_pos - dir * head + perp * head * 0.5f;
        Vector3f right = end_pos - dir * head - perp * head * 0.5f;

        glVertex3f(end_pos.x(), end_pos.y(), end_pos.z());
        glVertex3f(left.x(),    left.y(),    left.z());
        glVertex3f(end_pos.x(), end_pos.y(), end_pos.z());
        glVertex3f(right.x(),   right.y(),   right.z());
    }
    glEnd();

    glLineWidth(1.0f);
}

void InertiaDebugWidget::paintGL() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glMatrixMode(GL_PROJECTION);
    glLoadMatrixf(m_proj.constData());

    QMatrix4x4 view;
    view.lookAt(
        QVector3D(0.0f, 0.0f, 3.0f),
        QVector3D(0.0f, 0.0f, 0.0f),
        QVector3D(0.0f, 1.0f, 0.0f)
        );

    QMatrix4x4 model;
    model.rotate(10.0f, 1.0f, 0.0f, 0.0f);
    model.rotate(150.0f, 0.0f, 1.0f, 0.0f);

    QMatrix4x4 mv = view * model;

    glMatrixMode(GL_MODELVIEW);
    glLoadMatrixf(mv.constData());

    drawMesh(m_mesh);

    // green  = body rotational inertia (M_b top-left)
    drawPrincipalAxes(m_body_rot_axes,    m_body_rot_vals,    0.0f, 1.0f, 0.0f);

    // yellow = added rotational inertia (M_f top-left)
    drawPrincipalAxes(m_fluid_rot_axes,   m_fluid_rot_vals,   1.0f, 1.0f, 0.0f);

    // cyan   = added translational mass (M_f bottom-right)
    drawPrincipalAxes(m_fluid_trans_axes, m_fluid_trans_vals, 0.0f, 0.8f, 1.0f);
}

#ifdef __APPLE__
#pragma clang diagnostic pop
#endif
