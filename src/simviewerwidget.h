#pragma once

#ifdef __APPLE__
#define GL_SILENCE_DEPRECATION
#endif

#include <QOpenGLWidget>
#include <QElapsedTimer>
#include <QTimer>
#include <Eigen/Dense>
#include <vector>
#include <string>

class Camera;
class Shader;
class Shape;

class SimViewerWidget : public QOpenGLWidget
{
    Q_OBJECT

public:
    SimViewerWidget(std::vector<std::vector<Eigen::Vector3d>>& frames,
                    const std::vector<Eigen::Vector3i>& faces,
                    const std::vector<Eigen::Vector2f>& texcoords,
                    const std::vector<Eigen::Vector3i>& face_texcoord_ids,
                    const std::string& texturePath,
                    QWidget *parent = nullptr);
    ~SimViewerWidget();

private:
    void initializeGL() override;
    void paintGL() override;
    void resizeGL(int w, int h) override;

    void mousePressEvent(QMouseEvent *event) override;
    void mouseMoveEvent(QMouseEvent *event) override;
    void mouseReleaseEvent(QMouseEvent *event) override;
    void wheelEvent(QWheelEvent *event) override;
    void keyPressEvent(QKeyEvent *event) override;
    void keyReleaseEvent(QKeyEvent *event) override;

private slots:
    void tick();

private:
    std::vector<std::vector<Eigen::Vector3d>>& m_frames;
    std::vector<Eigen::Vector3i> m_faces;
    std::vector<Eigen::Vector2f> m_texcoords;
    std::vector<Eigen::Vector3i> m_face_texcoord_ids;
    std::string m_texturePath;

    int m_currentFrame = 0;
    bool m_paused = false;

    Camera *m_camera = nullptr;
    Shader *m_shader = nullptr;
    Shape *m_shape = nullptr;
    Shape *m_ground = nullptr;

    QTimer m_playbackTimer;

    bool m_capture = false;
    int m_lastX = 0;
    int m_lastY = 0;

    float m_forward = 0;
    float m_sideways = 0;
    float m_vertical = 0;

    void initGround();
};
