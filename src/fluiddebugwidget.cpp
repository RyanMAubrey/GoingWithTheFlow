#include "fluiddebugwidget.h"

#include <QDebug>
#include <QDir>
#include <QCoreApplication>
#include <QImage>

#include <algorithm>
#include <limits>
#include <filesystem>

#ifdef __APPLE__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif

using namespace Eigen;

FluidDebugWidget::FluidDebugWidget(QWidget *parent)
    : QOpenGLWidget(parent)
{
}

GLuint FluidDebugWidget::loadTextureQt(const std::string& path) {
    QImage img(QString::fromStdString(path));

    if (img.isNull()) {
        qDebug() << "Failed to load texture:" << QString::fromStdString(path);
        return 0;
    }

    img = img.convertToFormat(QImage::Format_RGBA8888).mirrored();

    GLuint tex = 0;
    glGenTextures(1, &tex);
    glBindTexture(GL_TEXTURE_2D, tex);

    glTexImage2D(
        GL_TEXTURE_2D,
        0,
        GL_RGBA,
        img.width(),
        img.height(),
        0,
        GL_RGBA,
        GL_UNSIGNED_BYTE,
        img.bits()
        );

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

    glBindTexture(GL_TEXTURE_2D, 0);

    qDebug() << "Loaded texture:" << QString::fromStdString(path)
             << "id:" << tex
             << "size:" << img.width() << img.height();

    return tex;
}

void FluidDebugWidget::loadMeshTextures(const TriMesh& mesh, const std::string& objPath) {
    std::filesystem::path objFile(objPath);
    std::string objDir = objFile.parent_path().string();

    for (const Material& mat : mesh.materials) {
        if (mat.map_kd.empty()) {
            continue;
        }

        if (m_textures.count(mat.name)) {
            continue;
        }

        std::string texPath;

        if (objDir.empty()) {
            texPath = mat.map_kd;
        } else {
            texPath = objDir + "/" + mat.map_kd;
        }

        GLuint tex = loadTextureQt(texPath);
        m_textures[mat.name] = tex;

        qDebug() << "Material:" << QString::fromStdString(mat.name)
                 << "map_Kd:" << QString::fromStdString(mat.map_kd);
    }
}

void FluidDebugWidget::loadFrames() {
    qDebug() << "Current working dir:" << QDir::currentPath();
    qDebug() << "App dir:" << QCoreApplication::applicationDirPath();

    std::string file_k  = "turtle_poses/frame_24.obj";
    std::string file_k1 = "turtle_poses/frame_32.obj";

    Integrator integrator;
    Momentum m;

    std::vector<Vector3f> vertices_k, vertices_k1;
    std::vector<Vector3i> faces_k, faces_k1;
    std::vector<Edge> edges_k, edges_k1;
    std::vector<float> face_areas_k, face_areas_k1;
    std::vector<Vector3f> face_normals_k, face_normals_k1;

    // One call now loads BOTH:
    // - full TriMesh for drawing/materials/textures
    // - physics arrays for momentum
    m_mesh_k = integrator.LoadPose(
        file_k,
        vertices_k,
        faces_k,
        edges_k,
        face_areas_k,
        face_normals_k
        );

    m_mesh_k1 = integrator.LoadPose(
        file_k1,
        vertices_k1,
        faces_k1,
        edges_k1,
        face_areas_k1,
        face_normals_k1
        );

    loadMeshTextures(m_mesh_k, file_k);

    qDebug() << "frame k verts:" << m_mesh_k.vertices.size()
             << "faces:" << m_mesh_k.faces.size()
             << "materials:" << m_mesh_k.materials.size()
             << "texcoords:" << m_mesh_k.texcoords.size()
             << "face tex ids:" << m_mesh_k.face_texcoord_ids.size()
             << "face material ids:" << m_mesh_k.face_material_ids.size();

    if (m_mesh_k.vertices.size() != m_mesh_k1.vertices.size()) {
        qDebug() << "ERROR: frames do not have matching vertex counts";
        return;
    }

    for (int i = 0; i < m_mesh_k.materials.size(); i++) {
        const Material& mat = m_mesh_k.materials[i];

        qDebug() << "material" << i
                 << QString::fromStdString(mat.name)
                 << "Kd:"
                 << mat.kd.x() << mat.kd.y() << mat.kd.z()
                 << "map_Kd:"
                 << QString::fromStdString(mat.map_kd);
    }

    float h = 7.0f;
    float rho_f = 1.0f;
    float delta = integrator.CalculateDelta(
        faces_k,
        edges_k,
        vertices_k,
        face_areas_k
        );

    qDebug() << "delta:" << delta;

    m_fluidMomentum = m.calc_fluid_momentum(
        vertices_k,
        vertices_k1,
        faces_k,
        face_areas_k,
        face_normals_k,
        edges_k,
        rho_f,
        h,
        delta
        );

    qDebug() << "Fluid momentum:"
             << m_fluidMomentum[0] << m_fluidMomentum[1] << m_fluidMomentum[2]
             << m_fluidMomentum[3] << m_fluidMomentum[4] << m_fluidMomentum[5];

    computeViewTransform();
}

void FluidDebugWidget::computeViewTransform() {
    if (m_mesh_k.vertices.empty()) {
        m_center = Vector3f::Zero();
        m_scale = 1.0f;
        return;
    }

    Vector3f minV(
        std::numeric_limits<float>::max(),
        std::numeric_limits<float>::max(),
        std::numeric_limits<float>::max()
        );

    Vector3f maxV(
        -std::numeric_limits<float>::max(),
        -std::numeric_limits<float>::max(),
        -std::numeric_limits<float>::max()
        );

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

Vector3f FluidDebugWidget::toViewSpace(const Vector3f& v) const {
    return (v - m_center) * m_scale;
}

void FluidDebugWidget::initializeGL() {
    initializeOpenGLFunctions();

    glEnable(GL_DEPTH_TEST);
    glDisable(GL_CULL_FACE);
    glDisable(GL_LIGHTING);

    glClearColor(0.1f, 0.1f, 0.1f, 1.0f);

    loadFrames();
}

void FluidDebugWidget::resizeGL(int w, int h) {
    glViewport(0, 0, w, h);

    m_proj.setToIdentity();
    m_proj.perspective(
        45.0f,
        float(w) / float(std::max(1, h)),
        0.1f,
        100.0f
        );
}

void FluidDebugWidget::drawMesh(const TriMesh& mesh) {
    for (int f = 0; f < mesh.faces.size(); f++) {
        int mat_id = -1;

        if (f < mesh.face_material_ids.size()) {
            mat_id = mesh.face_material_ids[f];
        }

        GLuint tex = 0;
        Vector3f kd(0.75f, 0.75f, 0.75f);

        if (mat_id >= 0 && mat_id < mesh.materials.size()) {
            const Material& mat = mesh.materials[mat_id];
            kd = mat.kd;

            if (!mat.map_kd.empty() && m_textures.count(mat.name)) {
                tex = m_textures[mat.name];
            }
        }

        bool hasTexture = tex != 0;

        if (hasTexture) {
            glEnable(GL_TEXTURE_2D);
            glBindTexture(GL_TEXTURE_2D, tex);

            // VERY IMPORTANT:
            // Do not use Kd here. Kd may be black.
            glColor3f(1.0f, 1.0f, 1.0f);
        } else {
            glDisable(GL_TEXTURE_2D);
            glBindTexture(GL_TEXTURE_2D, 0);
            glColor3f(kd.x(), kd.y(), kd.z());
        }

        const Vector3i& face = mesh.faces[f];

        Vector3i texFace(-1, -1, -1);
        if (f < mesh.face_texcoord_ids.size()) {
            texFace = mesh.face_texcoord_ids[f];
        }

        glBegin(GL_TRIANGLES);

        for (int corner = 0; corner < 3; corner++) {
            int vi = face[corner];
            int ti = texFace[corner];

            if (hasTexture && ti >= 0 && ti < mesh.texcoords.size()) {
                Vector2f uv = mesh.texcoords[ti];
                glTexCoord2f(uv.x(), uv.y());
            }

            Vector3f v = toViewSpace(mesh.vertices[vi]);
            glVertex3f(v.x(), v.y(), v.z());
        }

        glEnd();
    }

    glBindTexture(GL_TEXTURE_2D, 0);
    glDisable(GL_TEXTURE_2D);
}

void FluidDebugWidget::drawMotionLines(const TriMesh& mesh_k, const TriMesh& mesh_k1) {
    int n = std::min((int)mesh_k.vertices.size(), (int)mesh_k1.vertices.size());

    glDisable(GL_TEXTURE_2D);
    glColor3f(1.0f, 0.0f, 0.0f);
    glLineWidth(1.0f);

    glBegin(GL_LINES);

    for (int i = 0; i < n; i++) {
        Vector3f a = toViewSpace(mesh_k.vertices[i]);
        Vector3f b = toViewSpace(mesh_k1.vertices[i]);

        Vector3f dir = b - a;
        float len = dir.norm();

        if (len < 1e-6f) {
            continue;
        }

        dir.normalize();

        glVertex3f(a.x(), a.y(), a.z());
        glVertex3f(b.x(), b.y(), b.z());

        float arrowSize = 0.01f;
        Vector3f perp = dir.unitOrthogonal();

        Vector3f left  = b - dir * arrowSize + perp * arrowSize * 0.5f;
        Vector3f right = b - dir * arrowSize - perp * arrowSize * 0.5f;

        glVertex3f(b.x(), b.y(), b.z());
        glVertex3f(left.x(), left.y(), left.z());

        glVertex3f(b.x(), b.y(), b.z());
        glVertex3f(right.x(), right.y(), right.z());
    }

    glEnd();
}

void FluidDebugWidget::drawFluidMomentum() {
    glDisable(GL_TEXTURE_2D);

    Vector3f center = toViewSpace(m_center);

    Vector3f angular(
        m_fluidMomentum[0],
        m_fluidMomentum[1],
        m_fluidMomentum[2]
        );

    Vector3f linear(
        m_fluidMomentum[3],
        m_fluidMomentum[4],
        m_fluidMomentum[5]
        );

    float linearScale = 0.0005f;
    float angularScale = 0.00001f;

    auto drawArrow = [&](const Vector3f& start,
                         const Vector3f& vec,
                         float scale,
                         float r, float g, float b) {
        Vector3f end = start + vec * scale;

        Vector3f dir = end - start;
        float len = dir.norm();

        if (len < 1e-6f) {
            return;
        }

        dir.normalize();

        glColor3f(r, g, b);
        glLineWidth(4.0f);

        glBegin(GL_LINES);

        glVertex3f(start.x(), start.y(), start.z());
        glVertex3f(end.x(), end.y(), end.z());

        float arrowSize = 0.15f * len;
        arrowSize = std::clamp(arrowSize, 0.04f, 0.12f);

        Vector3f perp = dir.unitOrthogonal();

        Vector3f left  = end - dir * arrowSize + perp * arrowSize * 0.5f;
        Vector3f right = end - dir * arrowSize - perp * arrowSize * 0.5f;

        glVertex3f(end.x(), end.y(), end.z());
        glVertex3f(left.x(), left.y(), left.z());

        glVertex3f(end.x(), end.y(), end.z());
        glVertex3f(right.x(), right.y(), right.z());

        glEnd();
    };

    drawArrow(center, linear, linearScale, 0.0f, 0.5f, 1.0f);
    drawArrow(center, angular, angularScale, 0.0f, 1.0f, 0.0f);

    glLineWidth(1.0f);
}

void FluidDebugWidget::paintGL() {
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
    drawFluidMomentum();
}

#ifdef __APPLE__
#pragma clang diagnostic pop
#endif
