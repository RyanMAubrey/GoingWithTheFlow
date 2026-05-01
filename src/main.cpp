#include <QApplication>
#include <QSurfaceFormat>

#include "integrator.h"
#include "simviewerwidget.h"

#include <iostream>
#include <filesystem>

int main(int argc, char *argv[])
{
    // Run simulation first (no Qt needed for this part)
    std::filesystem::create_directories("output");

    Integrator integrator;
    integrator.Simulate();

    std::cout << "Simulation complete. " << integrator.output_frames.size()
              << " frames stored for playback." << std::endl;

    // Now open the viewer to play back the results
    QApplication app(argc, argv);

    QSurfaceFormat fmt;
    fmt.setVersion(3, 3);
    fmt.setProfile(QSurfaceFormat::CoreProfile);
    QSurfaceFormat::setDefaultFormat(fmt);

    // Get the face topology for rendering
    std::vector<Eigen::Vector3i> faces = integrator.getSharedFaces();

    SimViewerWidget viewer(integrator.output_frames, faces);
    viewer.resize(1200, 800);
    viewer.setWindowTitle("Going with the Flow — Turtle Swimming");
    viewer.show();

    return app.exec();
}
