#include "bodydebugwidget.h"

#include <QDebug>
#include <QDir>
#include <QCoreApplication>
#include <algorithm>
#include <limits>
#include <vector> // needed

#ifdef __APPLE__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif

using namespace Eigen;

BodyDebugWidget::BodyDebugWidget(QWidget *parent)
    : QOpenGLWidget(parent)
{
}

void BodyDebugWidget::loadFrames() {
    qDebug() << "Current working dir:" << QDir::currentPath();
    qDebug() << "App dir:" << QCoreApplication::applicationDirPath();

    m_mesh_k  = load_obj("turtle_poses/frame_1.obj");
    m_mesh_k1 = load_obj("turtle_poses/frame_8.obj");

    qDebug() << "frame 1 verts:" << m_mesh_k.vertices.size()
             << "faces:" << m_mesh_k.faces.size();
    qDebug() << "frame 8 verts:" << m_mesh_k1.vertices.size()
             << "faces:" << m_mesh_k1.faces.size();

    if (m_mesh_k.vertices.size() != m_mesh_k1.vertices.size()) {
        qDebug() << "ERROR: frame_1 and frame_8 do not have matching vertex counts";
    }

    Momentum m;
    std::vector<float> mass_density(m_mesh_k.vertices.size(), 1.0f);
    float h = 7.0f;

    m_bodyMomentum = m.calc_body_momentum(
        m_mesh_k.vertices,
        m_mesh_k1.vertices,
        mass_density,
        h
        );

    qDebug() << "Body momentum:"
             << m_bodyMomentum[0] << m_bodyMomentum[1] << m_bodyMomentum[2]
             << m_bodyMomentum[3] << m_bodyMomentum[4] << m_bodyMomentum[5];

    computeViewTransform();
}

void BodyDebugWidget::computeViewTransform() {
    if (m_mesh_k.vertices.empty()) {
        m_center = Vector3f::Zero();
        m_scale = 1.0f;
        return;
    }

    Vector3f minV( std::numeric_limits<float>::max(),
                  std::numeric_limits<float>::max(),
                  std::numeric_limits<float>::max() );

    Vector3f maxV( -std::numeric_limits<float>::max(),
                  -std::numeric_limits<float>::max(),
                  -std::numeric_limits<float>::max() );

    for (const auto& v : m_mesh_k.vertices) {
        minV = minV.cwiseMin(v);
        maxV = maxV.cwiseMax(v);
    }

    m_center = 0.5f * (minV + maxV);

    Vector3f extent = maxV - minV;
    float maxExtent = std::max(extent.x(), std::max(extent.y(), extent.z()));

    if (maxExtent > 0.0f) {
        m_scale = 1.8f / maxExtent;
    } else {
        m_scale = 1.0f;
    }

    qDebug() << "center:" << m_center.x() << m_center.y() << m_center.z();
    qDebug() << "scale:" << m_scale;
}

Vector3f BodyDebugWidget::toViewSpace(const Vector3f& v) const {
    return (v - m_center) * m_scale;
}

void BodyDebugWidget::initializeGL() {
    initializeOpenGLFunctions();

    glEnable(GL_DEPTH_TEST);
    glDisable(GL_CULL_FACE);
    glClearColor(0.1f, 0.1f, 0.1f, 1.0f);

    loadFrames();
}

void BodyDebugWidget::resizeGL(int w, int h) {
    glViewport(0, 0, w, h);

    m_proj.setToIdentity();
    m_proj.perspective(45.0f, float(w) / float(std::max(1, h)), 0.1f, 100.0f);
}

void BodyDebugWidget::drawMesh(const TriMesh& mesh) {
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

void BodyDebugWidget::drawMotionLines(const TriMesh& mesh_k, const TriMesh& mesh_k1) {
    int n = std::min((int)mesh_k.vertices.size(), (int)mesh_k1.vertices.size());

    glColor3f(1.0f, 0.0f, 0.0f);
    glLineWidth(1.0f);

    glBegin(GL_LINES);
    for (int i = 0; i < n; i++) {
        Eigen::Vector3f a = toViewSpace(mesh_k.vertices[i]);
        Eigen::Vector3f b = toViewSpace(mesh_k1.vertices[i]);

        Eigen::Vector3f dir = b - a;
        float len = dir.norm();
        if (len < 1e-6f) continue;

        dir.normalize();

        glVertex3f(a.x(), a.y(), a.z());
        glVertex3f(b.x(), b.y(), b.z());

        float arrowSize = 0.01f;
        Eigen::Vector3f perp = dir.unitOrthogonal();

        Eigen::Vector3f left  = b - dir * arrowSize + perp * arrowSize * 0.5f;
        Eigen::Vector3f right = b - dir * arrowSize - perp * arrowSize * 0.5f;

        glVertex3f(b.x(), b.y(), b.z());
        glVertex3f(left.x(), left.y(), left.z());

        glVertex3f(b.x(), b.y(), b.z());
        glVertex3f(right.x(), right.y(), right.z());
    }
    glEnd();
}

void BodyDebugWidget::drawBodyMomentum() {
    Vector3f center = toViewSpace(m_center);

    Vector3f angular(
        m_bodyMomentum[0],
        m_bodyMomentum[1],
        m_bodyMomentum[2]
        );

    Vector3f linear(
        m_bodyMomentum[3],
        m_bodyMomentum[4],
        m_bodyMomentum[5]
        );

    float linearScale = 0.0007f;
    float angularScale = 0.00002f;

    auto drawArrow = [&](const Vector3f& start,
                         const Vector3f& vec,
                         float scale,
                         float r, float g, float b) {
        Vector3f end = start + vec * scale;

        Vector3f dir = end - start;
        float len = dir.norm();
        if (len < 1e-6f) return;

        dir.normalize();

        glColor3f(r, g, b);
        glLineWidth(4.0f);

        glBegin(GL_LINES);

        // main line
        glVertex3f(start.x(), start.y(), start.z());
        glVertex3f(end.x(), end.y(), end.z());

        // bigger arrowhead
        float arrowSize = 0.15f * len;
        arrowSize = std::clamp(arrowSize, 0.04f, 0.12f);

        Vector3f perp = dir.unitOrthogonal();

        Vector3f left  = end - dir * arrowSize + perp * arrowSize * 0.5f;
        Vector3f right = end - dir * arrowSize - perp * arrowSize * 0.5f;

        // left wing
        glVertex3f(end.x(), end.y(), end.z());
        glVertex3f(left.x(), left.y(), left.z());

        // right wing
        glVertex3f(end.x(), end.y(), end.z());
        glVertex3f(right.x(), right.y(), right.z());

        glEnd();
    };

    // blue = linear momentum
    drawArrow(center, linear, linearScale, 0.0f, 0.5f, 1.0f);

    // green = angular momentum
    drawArrow(center, angular, angularScale, 0.0f, 1.0f, 0.0f);

    glLineWidth(1.0f);
}

void BodyDebugWidget::paintGL() {
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

    drawMesh(m_mesh_k);
    drawMotionLines(m_mesh_k, m_mesh_k1);
    drawBodyMomentum();
}

#ifdef __APPLE__
#pragma clang diagnostic pop
#endif
