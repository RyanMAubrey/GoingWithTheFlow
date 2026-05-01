#include <QApplication>
#include <QSurfaceFormat>

#include "integrator.h"
#include "simviewerwidget.h"

#include <iostream>
#include <filesystem>

int main(int argc, char *argv[])
{
    // Run simulation first
    std::filesystem::create_directories("output");

    Integrator integrator;
    integrator.Simulate();

    std::cout << "Simulation complete. " << integrator.output_frames.size()
              << " frames stored for playback." << std::endl;

    // Open the viewer
    QApplication app(argc, argv);

    QSurfaceFormat fmt;
    fmt.setVersion(3, 3);
    fmt.setProfile(QSurfaceFormat::CoreProfile);
    QSurfaceFormat::setDefaultFormat(fmt);

    std::vector<Eigen::Vector3i> faces = integrator.getSharedFaces();

    SimViewerWidget viewer(integrator.output_frames, faces,
                           integrator.texcoords,
                           integrator.face_texcoord_ids,
                           integrator.texture_path);
    viewer.resize(1200, 800);
    viewer.setWindowTitle("Going with the Flow");
    viewer.show();

    return app.exec();
}
