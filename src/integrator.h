#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <mesh_loader.h>

#include <vector>
#include <Eigen/Dense>

using namespace Eigen;
using Vector6f = Eigen::Matrix<float, 6, 1>;
using Matrix6f = Eigen::Matrix<float, 6, 6>;

// Put two Vector3f halves into a Vector6f.
static inline Vector6f make_vec6(const Vector3f& rot, const Vector3f& trans) {
    Vector6f v;
    v << rot, trans;
    return v;
}

struct Edge {
    int i;
    int j;
    float alpha; // [-pi, pi]
};

class Integrator
{
public:
    Integrator();

    void Simulate();

    // After Simulate() completes, these hold the output for visualization.
    // Vertices are in world space (Vector3d for compatibility with Shape class).
    std::vector<std::vector<Eigen::Vector3d>> output_frames;
    const std::vector<Vector3i>& getSharedFaces() const { return shared_faces; }

    TriMesh LoadPose(const std::string& filepath, std::vector<Vector3f>& vertices, std::vector<Vector3i>& faces, std::vector<Edge>& edges,
                  std::vector<float>& face_areas, std::vector<Vector3f>& face_normals);

    void CalculateFaceAttributes(std::vector<Vector3i>& faces, std::vector<Vector3f>& gamma,
                                 std::vector<float>& face_areas, std::vector<Vector3f>& face_normals);
    float CalculateDelta(std::vector<Vector3i>& faces, std::vector<Edge>& edges, std::vector<Vector3f>& gamma, std::vector<float>& face_areas);

    // Load all pose OBJs in order.
    void LoadAllPoses();

private:
    int total_frames = 40;
    float stroke_duration = 1.0f;
    float pair_duration = stroke_duration / (total_frames - 1);  // time per pose pair
    int substeps = 10;                                           // substeps per pose pair
    float h = pair_duration / substeps;                         // integrator timestep
    int num_strokes = 5;
    float rho_f = 998.0f;
    int output_every = 5;                                       // write an OBJ every N substeps

    // Pose storage from LoadAllPoses().
    std::vector<std::vector<Vector3f>> all_vertices;    // arrays of vertices
    std::vector<Vector3i> shared_faces;
    std::vector<Edge> shared_edges;
    std::vector<float> shared_face_areas_scratch;        
    std::vector<Vector3f> shared_face_normals_scratch;   

    bool HasEdge(Vector3i& face, int i, int j);
};

#endif // INTEGRATOR_H
