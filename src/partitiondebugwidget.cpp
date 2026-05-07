#include "partitiondebugwidget.h"

#include <QDebug>
#include <QDir>
#include <QCoreApplication>
#include <QVector4D>
#include <algorithm>
#include <array>
#include <cstdint>
#include <limits>
#include <queue>
#include <unordered_map>

#ifdef __APPLE__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif

using namespace Eigen;

namespace {

uint64_t make_edge_key(int a, int b)
{
    if (a > b) {
        std::swap(a, b);
    }

    return (static_cast<uint64_t>(static_cast<uint32_t>(a)) << 32) |
           static_cast<uint32_t>(b);
}

bool ray_triangle_intersection(const Vector3f& ray_origin,
                               const Vector3f& ray_dir,
                               const Vector3f& a,
                               const Vector3f& b,
                               const Vector3f& c,
                               float& t_out)
{
    const float eps = 1e-6f;
    const Vector3f ab = b - a;
    const Vector3f ac = c - a;
    const Vector3f p = ray_dir.cross(ac);
    const float det = ab.dot(p);

    if (std::abs(det) < eps) {
        return false;
    }

    const float inv_det = 1.0f / det;
    const Vector3f s = ray_origin - a;
    const float u = inv_det * s.dot(p);
    if (u < 0.0f || u > 1.0f) {
        return false;
    }

    const Vector3f q = s.cross(ab);
    const float v = inv_det * ray_dir.dot(q);
    if (v < 0.0f || u + v > 1.0f) {
        return false;
    }

    const float t = inv_det * ac.dot(q);
    if (t <= eps) {
        return false;
    }

    t_out = t;
    return true;
}

} // namespace

PartitionDebugWidget::PartitionDebugWidget(int frame, QWidget *parent)
    : QOpenGLWidget(parent), m_frame(frame)
{
    setFocusPolicy(Qt::StrongFocus);
    setMouseTracking(true);
}

// One distinct color per region. Index matches the TurtleRegion enum.
static const float REGION_COLORS[NUM_REGIONS][3] = {
    { 1.0f, 0.85f, 0.0f },  // REGION_HEAD               -> yellow
    { 0.2f, 1.0f,  0.2f },  // REGION_LEFT_FRONT_FLIPPER  -> green
    { 0.0f, 0.7f,  1.0f },  // REGION_RIGHT_FRONT_FLIPPER -> cyan
    { 1.0f, 0.5f,  0.0f },  // REGION_LEFT_REAR_FLIPPER   -> orange
    { 1.0f, 0.2f,  0.8f },  // REGION_RIGHT_REAR_FLIPPER  -> magenta
    { 0.72f, 0.72f, 0.72f }, // REGION_SHELL              -> light gray
};

void PartitionDebugWidget::loadFrames() {
    qDebug() << "Current working dir:" << QDir::currentPath();
    qDebug() << "App dir:" << QCoreApplication::applicationDirPath();

    const std::string path = "turtle_poses/frame_" + std::to_string(m_frame) + ".obj";
    m_mesh = load_obj(path);

    qDebug() << "frame" << m_frame
             << "verts:" << m_mesh.vertices.size()
             << "faces:" << m_mesh.faces.size();

    m_region = partition_mesh(m_mesh.vertices, m_mesh.faces);
    m_auto_region = m_region;
    buildFaceAdjacency();

    qDebug() << "Partition editor controls:";
    qDebug() << "  Keys 1-6: select region";
    qDebug() << "  Left drag: paint selected region";
    qDebug() << "  Right drag: paint shell";
    qDebug() << "  Shift + drag: orbit camera";
    qDebug() << "  Alt/Option + drag: pan camera";
    qDebug() << "  Two-finger scroll: zoom";
    qDebug() << "  [ and ]: shrink/grow brush steps";
    qDebug() << "  S: save labels to" << QString::fromStdString(partition_labels_path());
    qDebug() << "  L: reload saved labels";
    qDebug() << "  A: reset to labels loaded at startup";

    // Print how many faces ended up in each region — useful for tuning thresholds.
    int counts[NUM_REGIONS] = {0};
    for (int r : m_region) {
        if (r >= 0 && r < NUM_REGIONS) counts[r]++;
    }
    qDebug() << "Region face counts:";
    qDebug() << "  HEAD:               " << counts[REGION_HEAD];
    qDebug() << "  LEFT_FRONT_FLIPPER: " << counts[REGION_LEFT_FRONT_FLIPPER];
    qDebug() << "  RIGHT_FRONT_FLIPPER:" << counts[REGION_RIGHT_FRONT_FLIPPER];
    qDebug() << "  LEFT_REAR_FLIPPER:  " << counts[REGION_LEFT_REAR_FLIPPER];
    qDebug() << "  RIGHT_REAR_FLIPPER: " << counts[REGION_RIGHT_REAR_FLIPPER];
    qDebug() << "  SHELL:              " << counts[REGION_SHELL];

    computeViewTransform();
}

void PartitionDebugWidget::computeViewTransform() {
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
}

Vector3f PartitionDebugWidget::toViewSpace(const Vector3f& v) const {
    return (v - m_center) * m_scale;
}

void PartitionDebugWidget::buildFaceAdjacency()
{
    m_face_adjacency.assign(m_mesh.faces.size(), {});
    std::unordered_map<uint64_t, int> edge_owner;
    edge_owner.reserve(m_mesh.faces.size() * 3);

    for (int f = 0; f < static_cast<int>(m_mesh.faces.size()); f++) {
        const auto& face = m_mesh.faces[f];
        const std::array<std::pair<int, int>, 3> edges = {{
            {face[0], face[1]},
            {face[1], face[2]},
            {face[2], face[0]}
        }};

        for (const auto& edge : edges) {
            const uint64_t key = make_edge_key(edge.first, edge.second);
            auto it = edge_owner.find(key);
            if (it == edge_owner.end()) {
                edge_owner.emplace(key, f);
            } else {
                const int other = it->second;
                m_face_adjacency[f].push_back(other);
                m_face_adjacency[other].push_back(f);
            }
        }
    }
}

int PartitionDebugWidget::pickFace(const QPoint& pos) const
{
    if (m_mesh.faces.empty()) {
        return -1;
    }

    const float x = (2.0f * static_cast<float>(pos.x())) / std::max(1, width()) - 1.0f;
    const float y = 1.0f - (2.0f * static_cast<float>(pos.y())) / std::max(1, height());

    QMatrix4x4 inv_mvp = (m_proj * m_model_view).inverted();
    const QVector4D near_clip(x, y, -1.0f, 1.0f);
    const QVector4D far_clip(x, y, 1.0f, 1.0f);

    QVector4D near_obj4 = inv_mvp * near_clip;
    QVector4D far_obj4 = inv_mvp * far_clip;
    if (std::abs(near_obj4.w()) < 1e-6f || std::abs(far_obj4.w()) < 1e-6f) {
        return -1;
    }

    near_obj4 /= near_obj4.w();
    far_obj4 /= far_obj4.w();

    const Vector3f ray_origin(near_obj4.x(), near_obj4.y(), near_obj4.z());
    Vector3f ray_dir(far_obj4.x() - near_obj4.x(),
                     far_obj4.y() - near_obj4.y(),
                     far_obj4.z() - near_obj4.z());
    if (ray_dir.norm() < 1e-6f) {
        return -1;
    }
    ray_dir.normalize();

    int best_face = -1;
    float best_t = std::numeric_limits<float>::max();
    for (int f = 0; f < static_cast<int>(m_mesh.faces.size()); f++) {
        const auto& face = m_mesh.faces[f];
        const Vector3f a = toViewSpace(m_mesh.vertices[face[0]]);
        const Vector3f b = toViewSpace(m_mesh.vertices[face[1]]);
        const Vector3f c = toViewSpace(m_mesh.vertices[face[2]]);

        float t = 0.0f;
        if (ray_triangle_intersection(ray_origin, ray_dir, a, b, c, t) && t < best_t) {
            best_t = t;
            best_face = f;
        }
    }

    return best_face;
}

void PartitionDebugWidget::paintFromFace(int face_index, int region_id)
{
    if (face_index < 0 || face_index >= static_cast<int>(m_region.size())) {
        return;
    }

    std::queue<std::pair<int, int>> q;
    std::vector<char> visited(m_region.size(), false);
    q.push({face_index, 0});
    visited[face_index] = true;

    while (!q.empty()) {
        const auto [face, depth] = q.front();
        q.pop();
        m_region[face] = region_id;

        if (depth >= m_brush_steps) {
            continue;
        }

        for (int neighbor : m_face_adjacency[face]) {
            if (!visited[neighbor]) {
                visited[neighbor] = true;
                q.push({neighbor, depth + 1});
            }
        }
    }

    update();
}

void PartitionDebugWidget::saveLabels() const
{
    if (save_partition_labels(partition_labels_path(), m_region)) {
        qDebug() << "Saved labels to" << QString::fromStdString(partition_labels_path());
    } else {
        qDebug() << "Failed to save labels to" << QString::fromStdString(partition_labels_path());
    }
}

void PartitionDebugWidget::reloadSavedLabels()
{
    const std::vector<int> loaded = load_partition_labels(partition_labels_path(),
                                                          static_cast<int>(m_mesh.faces.size()));
    if (!loaded.empty()) {
        m_region = loaded;
        qDebug() << "Reloaded saved labels from" << QString::fromStdString(partition_labels_path());
    } else {
        qDebug() << "No saved labels found at" << QString::fromStdString(partition_labels_path());
    }
    update();
}

void PartitionDebugWidget::initializeGL() {
    initializeOpenGLFunctions();

    glEnable(GL_DEPTH_TEST);
    glDisable(GL_CULL_FACE);
    glClearColor(0.1f, 0.1f, 0.1f, 1.0f);

    loadFrames();
}

void PartitionDebugWidget::resizeGL(int w, int h) {
    glViewport(0, 0, w, h);

    m_proj.setToIdentity();
    m_proj.perspective(45.0f, float(w) / float(std::max(1, h)), 0.1f, 100.0f);
}

void PartitionDebugWidget::drawColoredMesh() {
    glBegin(GL_TRIANGLES);
    for (int f = 0; f < (int)m_mesh.faces.size(); f++) {
        int r = (f < (int)m_region.size()) ? m_region[f] : REGION_SHELL;
        if (r < 0 || r >= NUM_REGIONS) r = REGION_SHELL;

        const float* col = REGION_COLORS[r];
        glColor3f(col[0], col[1], col[2]);

        const auto& face = m_mesh.faces[f];
        const auto a = toViewSpace(m_mesh.vertices[face[0]]);
        const auto b = toViewSpace(m_mesh.vertices[face[1]]);
        const auto c = toViewSpace(m_mesh.vertices[face[2]]);

        glVertex3f(a.x(), a.y(), a.z());
        glVertex3f(b.x(), b.y(), b.z());
        glVertex3f(c.x(), c.y(), c.z());
    }
    glEnd();
}

void PartitionDebugWidget::paintGL() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glMatrixMode(GL_PROJECTION);
    glLoadMatrixf(m_proj.constData());

    QMatrix4x4 view;
    view.lookAt(
        QVector3D(0.0f, 0.0f, m_camera_distance),
        QVector3D(0.0f, 0.0f, 0.0f),
        QVector3D(0.0f, 1.0f, 0.0f)
        );

    QMatrix4x4 model;
    model.translate(m_pan);
    model.rotate(m_pitch_degrees, 1.0f, 0.0f, 0.0f);
    model.rotate(m_yaw_degrees, 0.0f, 1.0f, 0.0f);

    QMatrix4x4 mv = view * model;
    m_model_view = mv;

    glMatrixMode(GL_MODELVIEW);
    glLoadMatrixf(mv.constData());

    drawColoredMesh();
}

void PartitionDebugWidget::mousePressEvent(QMouseEvent *event)
{
    m_last_mouse_pos = event->pos();

    if ((event->modifiers() & Qt::ShiftModifier) != 0) {
        m_is_rotating = true;
        return;
    }
    if ((event->modifiers() & Qt::AltModifier) != 0) {
        m_is_panning = true;
        return;
    }

    const int face = pickFace(event->pos());
    if (face < 0) {
        return;
    }

    const int region = (event->button() == Qt::RightButton) ? REGION_SHELL : m_active_region;
    m_is_painting = true;
    paintFromFace(face, region);
}

void PartitionDebugWidget::mouseMoveEvent(QMouseEvent *event)
{
    const QPoint delta = event->pos() - m_last_mouse_pos;
    m_last_mouse_pos = event->pos();

    if (m_is_rotating) {
        m_yaw_degrees += 0.4f * static_cast<float>(delta.x());
        m_pitch_degrees += 0.4f * static_cast<float>(delta.y());
        m_pitch_degrees = std::clamp(m_pitch_degrees, -89.0f, 89.0f);
        update();
        return;
    }

    if (m_is_panning) {
        const float scale = 0.0035f * std::max(1.0f, m_camera_distance);
        m_pan += QVector3D(scale * static_cast<float>(delta.x()),
                           -scale * static_cast<float>(delta.y()),
                           0.0f);
        update();
        return;
    }

    if (!(event->buttons() & (Qt::LeftButton | Qt::RightButton))) {
        return;
    }

    const int face = pickFace(event->pos());
    if (face < 0) {
        return;
    }

    const int region = (event->buttons() & Qt::RightButton) ? REGION_SHELL : m_active_region;
    paintFromFace(face, region);
}

void PartitionDebugWidget::mouseReleaseEvent(QMouseEvent *event)
{
    Q_UNUSED(event);
    m_is_painting = false;
    m_is_rotating = false;
    m_is_panning = false;
}

void PartitionDebugWidget::wheelEvent(QWheelEvent *event)
{
    const QPoint num_degrees = event->angleDelta();
    if (num_degrees.isNull()) {
        event->ignore();
        return;
    }

    const float zoom_step = -0.0015f * static_cast<float>(num_degrees.y());
    m_camera_distance = std::clamp(m_camera_distance + zoom_step, 1.2f, 12.0f);
    update();
    event->accept();
}

void PartitionDebugWidget::keyPressEvent(QKeyEvent *event)
{
    switch (event->key()) {
    case Qt::Key_1: m_active_region = REGION_HEAD; break;
    case Qt::Key_2: m_active_region = REGION_LEFT_FRONT_FLIPPER; break;
    case Qt::Key_3: m_active_region = REGION_RIGHT_FRONT_FLIPPER; break;
    case Qt::Key_4: m_active_region = REGION_LEFT_REAR_FLIPPER; break;
    case Qt::Key_5: m_active_region = REGION_RIGHT_REAR_FLIPPER; break;
    case Qt::Key_6: m_active_region = REGION_SHELL; break;
    case Qt::Key_BracketLeft:
        m_brush_steps = std::max(0, m_brush_steps - 1);
        qDebug() << "Brush steps:" << m_brush_steps;
        break;
    case Qt::Key_BracketRight:
        m_brush_steps = std::min(10, m_brush_steps + 1);
        qDebug() << "Brush steps:" << m_brush_steps;
        break;
    case Qt::Key_S:
        saveLabels();
        break;
    case Qt::Key_L:
        reloadSavedLabels();
        break;
    case Qt::Key_A:
        m_region = m_auto_region;
        qDebug() << "Reset labels to startup state";
        update();
        break;
    default:
        QOpenGLWidget::keyPressEvent(event);
        return;
    }

    qDebug() << "Active region:" << m_active_region;
}

#ifdef __APPLE__
#pragma clang diagnostic pop
#endif
