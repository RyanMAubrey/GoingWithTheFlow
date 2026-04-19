#pragma once

#include <QOpenGLWidget>
#include <QOpenGLFunctions>
#include <QMatrix4x4>
#include <Eigen/Dense>

#include "mesh_loader.h"
#include "momentum.h"

class BodyDebugWidget : public QOpenGLWidget, protected QOpenGLFunctions
{
    Q_OBJECT

public:
    explicit BodyDebugWidget(QWidget *parent = nullptr);

protected:
    void initializeGL() override;
    void resizeGL(int w, int h) override;
    void paintGL() override;

private:
    TriMesh m_mesh_k;
    TriMesh m_mesh_k1;
    QMatrix4x4 m_proj;

    Eigen::Vector3f m_center = Eigen::Vector3f::Zero();
    float m_scale = 1.0f;
    Vector6f m_bodyMomentum = Vector6f::Zero();

    void loadFrames();
    void computeViewTransform();
    Eigen::Vector3f toViewSpace(const Eigen::Vector3f& v) const;

    void drawMesh(const TriMesh& mesh);
    void drawMotionLines(const TriMesh& mesh_k, const TriMesh& mesh_k1);
    void drawBodyMomentum();
};
