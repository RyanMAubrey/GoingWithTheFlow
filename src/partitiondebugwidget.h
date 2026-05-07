#pragma once

#include <QOpenGLWidget>
#include <QOpenGLFunctions>
#include <QMatrix4x4>
#include <QMouseEvent>
#include <QKeyEvent>
#include <QWheelEvent>
#include <Eigen/Dense>

#include "mesh_loader.h"
#include "partition.h"

// Visualizes the output of partition_mesh() by coloring each triangle
// of the turtle according to its region label (head, four flippers, shell).
//
// Use this to sanity-check the AABB-based classification rules: if a
// flipper shows up gray (shell-colored), the rule thresholds need tuning.
class PartitionDebugWidget : public QOpenGLWidget, protected QOpenGLFunctions
{
    Q_OBJECT

public:
    explicit PartitionDebugWidget(int frame = 1, QWidget *parent = nullptr);

protected:
    void initializeGL() override;
    void resizeGL(int w, int h) override;
    void paintGL() override;
    void mousePressEvent(QMouseEvent *event) override;
    void mouseMoveEvent(QMouseEvent *event) override;
    void mouseReleaseEvent(QMouseEvent *event) override;
    void wheelEvent(QWheelEvent *event) override;
    void keyPressEvent(QKeyEvent *event) override;

private:
    int m_frame = 1;
    TriMesh m_mesh;
    std::vector<int> m_region;
    std::vector<int> m_auto_region;
    std::vector<std::vector<int>> m_face_adjacency;
    QMatrix4x4 m_proj;
    QMatrix4x4 m_model_view;

    Eigen::Vector3f m_center = Eigen::Vector3f::Zero();
    float m_scale = 1.0f;
    float m_yaw_degrees = 150.0f;
    float m_pitch_degrees = 10.0f;
    float m_camera_distance = 3.0f;
    QVector3D m_pan = QVector3D(0.0f, 0.0f, 0.0f);
    int m_active_region = REGION_HEAD;
    int m_brush_steps = 2;
    QPoint m_last_mouse_pos;
    bool m_is_painting = false;
    bool m_is_rotating = false;
    bool m_is_panning = false;

    void loadFrames();
    void computeViewTransform();
    Eigen::Vector3f toViewSpace(const Eigen::Vector3f& v) const;
    void buildFaceAdjacency();
    int pickFace(const QPoint& pos) const;
    void paintFromFace(int face_index, int region_id);
    void saveLabels() const;
    void reloadSavedLabels();

    void drawColoredMesh();
};
