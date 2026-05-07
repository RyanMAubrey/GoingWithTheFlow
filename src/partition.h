#ifndef PARTITION_H
#define PARTITION_H

#include <string>
#include <vector>
#include <Eigen/Dense>

// Region ID's for turle mesh
enum TurtleRegion {
    REGION_HEAD = 0,
    REGION_LEFT_FRONT_FLIPPER = 1,
    REGION_RIGHT_FRONT_FLIPPER = 2 ,
    REGION_LEFT_REAR_FLIPPER = 3 ,
    REGION_RIGHT_REAR_FLIPPER = 4,
    REGION_SHELL = 5,
    NUM_REGIONS = 6
};

std::vector<int> partition_mesh(const std::vector<Eigen::Vector3f>& verticies,
                                const std::vector<Eigen::Vector3i>& faces);

std::string partition_labels_path();

std::vector<int> load_partition_labels(const std::string& path, int expected_faces);
bool save_partition_labels(const std::string& path, const std::vector<int>& region);

#endif
