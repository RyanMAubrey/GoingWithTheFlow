#include <GL/glew.h>  // Must be first — before any Qt/GL headers
#include "simviewerwidget.h"
#include "graphics/camera.h"
#include "graphics/shader.h"
#include "graphics/shape.h"

#include <QApplication>
#include <QKeyEvent>
#include <iostream>

#define SPEED 1.5
#define ROTATE_SPEED 0.0025

using namespace Eigen;

SimViewerWidget::SimViewerWidget(std::vector<std::vector<Vector3d>>& frames,
                                 const std::vector<Vector3i>& faces,
                                 QWidget *parent)
    : QOpenGLWidget(parent),
      m_frames(frames),
      m_faces(faces)
{
    setMouseTracking(true);
    QApplication::setOverrideCursor(Qt::ArrowCursor);
    setFocusPolicy(Qt::StrongFocus);

    connect(&m_playbackTimer, SIGNAL(timeout()), this, SLOT(tick()));
}

SimViewerWidget::~SimViewerWidget()
{
    delete m_shader;
    delete m_camera;
    delete m_shape;
    delete m_ground;
}

void SimViewerWidget::initializeGL()
{
    glewExperimental = GL_TRUE;
    GLenum err = glewInit();
    if (err != GLEW_OK)
        fprintf(stderr, "GLEW init error: %s\n", glewGetErrorString(err));

    glClearColor(0.05f, 0.05f, 0.1f, 1.0f);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);

    m_shader = new Shader(":/resources/shaders/shader.vert",
                          ":/resources/shaders/shader.frag");

    m_camera = new Camera();
    m_shape = new Shape();
    m_ground = new Shape();

    // Initialize the shape with the first frame
    if (!m_frames.empty()) {
        m_shape->init(m_frames[0], m_faces);
        m_shape->setColor(0.4f, 0.7f, 0.45f);  // greenish turtle color
    }

    initGround();
    m_ground->setColor(0.2f, 0.3f, 0.6f);  // blue-ish water floor

    // Set up camera — turtle is ~0.6m long, swims in roughly +Y/-Z
    Vector3f eye(0.5f, 0.3f, 0.5f);
    Vector3f target(0.0f, 0.05f, 0.0f);
    m_camera->lookAt(eye, target);
    m_camera->setOrbitPoint(target);
    m_camera->toggleIsOrbiting();
    m_camera->setPerspective(45, width() / static_cast<float>(height()), 0.001f, 200.0f);

    m_playbackTimer.start(1000 / 30);  // 30 fps playback
}

void SimViewerWidget::paintGL()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    m_shader->bind();
    m_shader->setUniform("proj", m_camera->getProjection());
    m_shader->setUniform("view", m_camera->getView());

    m_shape->draw(m_shader);
    m_ground->draw(m_shader);

    m_shader->unbind();
}

void SimViewerWidget::resizeGL(int w, int h)
{
    glViewport(0, 0, w, h);
    m_camera->setAspect(static_cast<float>(w) / h);
}

// ================== Input handling ==================

void SimViewerWidget::mousePressEvent(QMouseEvent *event)
{
    m_capture = true;
    m_lastX = event->position().x();
    m_lastY = event->position().y();
}

void SimViewerWidget::mouseMoveEvent(QMouseEvent *event)
{
    if (!m_capture) return;

    int currX = event->position().x();
    int currY = event->position().y();
    int deltaX = currX - m_lastX;
    int deltaY = currY - m_lastY;

    if (deltaX == 0 && deltaY == 0) return;

    m_camera->rotate(deltaY * ROTATE_SPEED, -deltaX * ROTATE_SPEED);

    m_lastX = currX;
    m_lastY = currY;
}

void SimViewerWidget::mouseReleaseEvent(QMouseEvent *event)
{
    m_capture = false;
}

void SimViewerWidget::wheelEvent(QWheelEvent *event)
{
    float zoom = 1 - event->pixelDelta().y() * 0.1f / 120.f;
    m_camera->zoom(zoom);
}

void SimViewerWidget::keyPressEvent(QKeyEvent *event)
{
    if (event->isAutoRepeat()) return;

    switch (event->key())
    {
    case Qt::Key_W: m_forward  += SPEED; break;
    case Qt::Key_S: m_forward  -= SPEED; break;
    case Qt::Key_A: m_sideways -= SPEED; break;
    case Qt::Key_D: m_sideways += SPEED; break;
    case Qt::Key_F: m_vertical -= SPEED; break;
    case Qt::Key_R: m_vertical += SPEED; break;
    case Qt::Key_C: m_camera->toggleIsOrbiting(); break;
    case Qt::Key_T: m_shape->toggleWireframe(); break;
    case Qt::Key_Space: m_paused = !m_paused; break;
    case Qt::Key_Escape: QApplication::quit(); break;
    }
}

void SimViewerWidget::keyReleaseEvent(QKeyEvent *event)
{
    if (event->isAutoRepeat()) return;

    switch (event->key())
    {
    case Qt::Key_W: m_forward  -= SPEED; break;
    case Qt::Key_S: m_forward  += SPEED; break;
    case Qt::Key_A: m_sideways += SPEED; break;
    case Qt::Key_D: m_sideways -= SPEED; break;
    case Qt::Key_F: m_vertical += SPEED; break;
    case Qt::Key_R: m_vertical -= SPEED; break;
    }
}

// ================== Playback tick ==================

void SimViewerWidget::tick()
{
    // Move camera
    auto look = m_camera->getLook();
    look.y() = 0;
    look.normalize();
    Vector3f perp(-look.z(), 0, look.x());
    Vector3f moveVec = m_forward * look.normalized()
                     + m_sideways * perp.normalized()
                     + m_vertical * Vector3f::UnitY();
    float dt = 1.0f / 30.0f;
    moveVec *= dt;
    m_camera->move(moveVec);

    // Advance frame
    if (!m_paused && !m_frames.empty()) {
        m_shape->setVertices(m_frames[m_currentFrame]);
        m_currentFrame++;
        if (m_currentFrame >= (int)m_frames.size()) {
            m_currentFrame = 0;  // loop
        }
    }

    update();  // trigger paintGL
}

// ================== Ground plane ==================

void SimViewerWidget::initGround()
{
    std::vector<Vector3d> groundVerts;
    std::vector<Vector3i> groundFaces;
    groundVerts.emplace_back(-50, -0.5, -50);
    groundVerts.emplace_back(-50, -0.5, 50);
    groundVerts.emplace_back(50, -0.5, 50);
    groundVerts.emplace_back(50, -0.5, -50);
    groundFaces.emplace_back(0, 1, 2);
    groundFaces.emplace_back(0, 2, 3);
    m_ground->init(groundVerts, groundFaces);
}
