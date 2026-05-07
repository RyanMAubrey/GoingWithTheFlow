#include <QApplication>
#include <QSurfaceFormat>

#include "integrator.h"
#include "partitiondebugwidget.h"
#include "simviewerwidget.h"

#include <algorithm>
#include <filesystem>
#include <iostream>
#include <string>

namespace {

struct ComparisonMetrics {
    double mean_trajectory_diff = 0.0;
    double max_trajectory_diff = 0.0;
    double final_trajectory_diff = 0.0;
    Eigen::Vector3d global_displacement = Eigen::Vector3d::Zero();
    Eigen::Vector3d piecewise_displacement = Eigen::Vector3d::Zero();
};

ComparisonMetrics compare_centroid_trajectories(const Integrator& global_integrator,
                                                const Integrator& piecewise_integrator)
{
    ComparisonMetrics metrics;

    const auto& global = global_integrator.output_centroids;
    const auto& piecewise = piecewise_integrator.output_centroids;
    const int count = std::min(global.size(), piecewise.size());
    if (count == 0) {
        return metrics;
    }

    double sum = 0.0;
    for (int i = 0; i < count; i++) {
        const double d = (piecewise[i] - global[i]).norm();
        sum += d;
        metrics.max_trajectory_diff = std::max(metrics.max_trajectory_diff, d);
    }

    metrics.mean_trajectory_diff = sum / static_cast<double>(count);
    metrics.final_trajectory_diff = (piecewise[count - 1] - global[count - 1]).norm();
    metrics.global_displacement = global[count - 1] - global[0];
    metrics.piecewise_displacement = piecewise[count - 1] - piecewise[0];
    return metrics;
}

void print_metrics(const ComparisonMetrics& metrics)
{
    std::cout << "\n=== Global vs Piecewise Added Mass ===" << std::endl;
    std::cout << "Mean centroid trajectory diff: " << metrics.mean_trajectory_diff << std::endl;
    std::cout << "Max centroid trajectory diff:  " << metrics.max_trajectory_diff << std::endl;
    std::cout << "Final centroid diff:           " << metrics.final_trajectory_diff << std::endl;
    std::cout << "Global displacement:           "
              << metrics.global_displacement.transpose() << std::endl;
    std::cout << "Piecewise displacement:        "
              << metrics.piecewise_displacement.transpose() << std::endl;
}

void configure_surface_format(bool partition_mode)
{
    QSurfaceFormat fmt;
    if (partition_mode) {
        fmt.setVersion(2, 1);
        fmt.setProfile(QSurfaceFormat::NoProfile);
    } else {
        fmt.setVersion(3, 3);
        fmt.setProfile(QSurfaceFormat::CoreProfile);
    }
    QSurfaceFormat::setDefaultFormat(fmt);
}

} // namespace

int main(int argc, char *argv[])
{
    const std::string mode = (argc > 1) ? std::string(argv[1]) : "";
    std::filesystem::create_directories("output");
    const bool partition_mode = (mode == "partition");

    if (partition_mode) {
        configure_surface_format(true);
        QApplication app(argc, argv);
        PartitionDebugWidget viewer(20);
        viewer.resize(1000, 800);
        viewer.setWindowTitle("Partition Editor - Frame 20");
        viewer.show();
        return app.exec();
    }

    if (mode == "compare") {
        Integrator global_integrator;
        global_integrator.Simulate(AddedMassMode::Global);

        Integrator piecewise_integrator;
        piecewise_integrator.Simulate(AddedMassMode::Piecewise);

        std::cout << "Global simulation complete. " << global_integrator.output_frames.size()
                  << " frames stored for playback." << std::endl;
        std::cout << "Piecewise simulation complete. " << piecewise_integrator.output_frames.size()
                  << " frames stored for playback." << std::endl;
        print_metrics(compare_centroid_trajectories(global_integrator, piecewise_integrator));

        configure_surface_format(false);
        QApplication app(argc, argv);
        std::vector<Eigen::Vector3i> global_faces = global_integrator.getSharedFaces();
        std::vector<Eigen::Vector3i> piecewise_faces = piecewise_integrator.getSharedFaces();

        SimViewerWidget global_viewer(global_integrator.output_frames, global_faces,
                                      global_integrator.texcoords,
                                      global_integrator.face_texcoord_ids,
                                      global_integrator.texture_path);
        global_viewer.resize(900, 700);
        global_viewer.setWindowTitle("Going with the Flow - Global Added Mass");
        global_viewer.move(80, 80);
        global_viewer.show();

        SimViewerWidget piecewise_viewer(piecewise_integrator.output_frames, piecewise_faces,
                                         piecewise_integrator.texcoords,
                                         piecewise_integrator.face_texcoord_ids,
                                         piecewise_integrator.texture_path);
        piecewise_viewer.resize(900, 700);
        piecewise_viewer.setWindowTitle("Going with the Flow - Piecewise Added Mass");
        piecewise_viewer.move(1020, 80);
        piecewise_viewer.show();
        return app.exec();
    }

    const AddedMassMode sim_mode = (mode == "piecewise")
        ? AddedMassMode::Piecewise
        : AddedMassMode::Global;

    Integrator integrator;
    integrator.Simulate(sim_mode);

    std::cout << "Simulation complete. " << integrator.output_frames.size()
              << " frames stored for playback." << std::endl;

    configure_surface_format(false);
    QApplication app(argc, argv);
    std::vector<Eigen::Vector3i> faces = integrator.getSharedFaces();

    SimViewerWidget viewer(integrator.output_frames, faces,
                           integrator.texcoords,
                           integrator.face_texcoord_ids,
                           integrator.texture_path);
    viewer.resize(1200, 800);
    viewer.setWindowTitle(mode == "piecewise"
                              ? "Going with the Flow - Piecewise Added Mass"
                              : "Going with the Flow");
    viewer.show();

    return app.exec();
}
