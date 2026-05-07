#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <mesh_loader.h>
#include <partition.h>

#include <array>
#include <utility>
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

enum class AddedMassMode {
    Global,
    Piecewise
};

class Integrator
{
public:
    Integrator();

    void Simulate();
    void Simulate(AddedMassMode mode);

    // These hold the output for visualization
    std::vector<std::vector<Eigen::Vector3d>> output_frames;
    std::vector<Eigen::Vector3d> output_centroids;
    const std::vector<Vector3i>& getSharedFaces() const { return shared_faces; }

    // Texture data from the first frame (for rendering)
    std::vector<Eigen::Vector2f> texcoords;
    std::vector<Eigen::Vector3i> face_texcoord_ids;
    std::string texture_path;

    TriMesh LoadPose(const std::string& filepath, std::vector<Vector3f>& vertices, 
                        std::vector<Vector3i>& faces, std::vector<Edge>& edges,
                        std::vector<float>& face_areas, std::vector<Vector3f>& face_normals);

    void CalculateFaceAttributes(std::vector<Vector3i>& faces, std::vector<Vector3f>& gamma,
                                 std::vector<float>& face_areas, std::vector<Vector3f>& face_normals);
    
    float CalculateDelta(std::vector<Vector3i>& faces, std::vector<Edge>& edges, 
                        std::vector<Vector3f>& gamma, std::vector<float>& face_areas);

    // Load all pose OBJs in order.
    void LoadAllPoses();

private:

    // General scene parameters
    int total_frames = 40;
    float stroke_duration = 1.0f;
    float pair_duration = stroke_duration / (total_frames - 1);     // time per pose pair
    int substeps = 10;                                              // substeps per pose pair
    float h = pair_duration / substeps;                             // integrator timestep
    int num_strokes = 10;                                           // number of swimming strokes to simulate
    float rho_f = 998.0f;                                           // density of fluid (kg/m^3)
    int output_every = 5;                                           // save output every N frames

    // Pose storage from LoadAllPoses().
    std::vector<std::vector<Vector3f>> all_vertices;    
    std::vector<Vector3i> shared_faces;
    std::vector<int> shared_face_regions;
    std::vector<Edge> shared_edges;
    std::vector<std::pair<int, int>> shared_edge_faces;
    std::vector<float> shared_face_areas_scratch;        
    std::vector<Vector3f> shared_face_normals_scratch;   

    bool HasEdge(const Vector3i& face, int i, int j) const;
    void BuildEdgeFaceAdjacency();
    std::array<float, NUM_REGIONS> CalculateRegionDeltas(const std::vector<Vector3f>& gamma,
                                                         const std::vector<float>& face_areas,
                                                         float fallback_delta) const;
    std::vector<float> BuildFaceDeltas(const std::array<float, NUM_REGIONS>& region_deltas,
                                       float fallback_delta) const;

    // Turtle scene parameters
    // std::string animal = "turtle";
    // float total_mass = 0.25f;                                   // slight effect of gravity
    // Vector3f gravity = Vector3f(0.0f, -9.81f, 0.0f);            // standard gravity
    // Vector3f bg_flow = Vector3f::Zero();                        // no background flow
    // float angular_damping = 0.95f;                              // 1.0 = no damping, 0.9 = strong damping
    // float drag_scale = 0.2f;                                    // 1.0 = full flat-plate drag, lower = more streamlined

    // Worm scene parameters
    std::string animal = "worm";
    float total_mass = 0.0f;                                   // slight effect of gravity
    Vector3f gravity = Vector3f(0.0f, -9.81f, 0.0f);            // standard gravity
    Vector3f bg_flow = Vector3f::Zero();                        // no background flow
    float angular_damping = 1.0f;                              // 1.0 = no damping, 0.9 = strong damping
    float drag_scale = 1.0f;                                    // 1.0 = full flat-plate drag, lower = more streamlined
};

#endif // INTEGRATOR_H
