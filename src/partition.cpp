#include "partition.h"
#include "mesh_loader.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <limits>
#include <numeric>
#include <queue>
#include <unordered_map>

using namespace Eigen;

namespace {

struct FaceComponent {
    std::vector<int> faces;
    Vector3f centroid = Vector3f::Zero();
    float max_distance = 0.0f;
    float mean_distance = 0.0f;
    float score = 0.0f;
};

uint64_t make_edge_key(int a, int b)
{
    if (a > b) {
        std::swap(a, b);
    }

    return (static_cast<uint64_t>(static_cast<uint32_t>(a)) << 32) |
           static_cast<uint32_t>(b);
}

Vector3f compute_mesh_centroid(const std::vector<Vector3f>& vertices)
{
    Vector3f centroid = Vector3f::Zero();
    for (const auto& v : vertices) {
        centroid += v;
    }

    return centroid / std::max(1.0f, static_cast<float>(vertices.size()));
}

std::vector<Vector3f> compute_face_centroids(const std::vector<Vector3f>& vertices,
                                             const std::vector<Vector3i>& faces)
{
    std::vector<Vector3f> centroids(faces.size(), Vector3f::Zero());
    for (int f = 0; f < static_cast<int>(faces.size()); f++) {
        centroids[f] = (vertices[faces[f][0]] +
                        vertices[faces[f][1]] +
                        vertices[faces[f][2]]) / 3.0f;
    }

    return centroids;
}

std::vector<std::vector<int>> build_face_adjacency(const std::vector<Vector3i>& faces)
{
    std::vector<std::vector<int>> adjacency(faces.size());
    std::unordered_map<uint64_t, int> edge_owner;
    edge_owner.reserve(faces.size() * 3);

    for (int f = 0; f < static_cast<int>(faces.size()); f++) {
        const Vector3i& face = faces[f];
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
                adjacency[f].push_back(other);
                adjacency[other].push_back(f);
            }
        }
    }

    return adjacency;
}

std::vector<float> compute_face_distances(const std::vector<Vector3f>& centroids,
                                          const std::vector<std::vector<int>>& adjacency,
                                          int start_face)
{
    const int num_faces = static_cast<int>(centroids.size());
    const float inf = std::numeric_limits<float>::infinity();
    std::vector<float> distance(num_faces, inf);

    using QueueEntry = std::pair<float, int>;
    std::priority_queue<QueueEntry, std::vector<QueueEntry>, std::greater<QueueEntry>> pq;

    distance[start_face] = 0.0f;
    pq.push({0.0f, start_face});

    while (!pq.empty()) {
        const auto [curr_dist, face] = pq.top();
        pq.pop();

        if (curr_dist > distance[face]) {
            continue;
        }

        for (int neighbor : adjacency[face]) {
            const float edge_length = (centroids[face] - centroids[neighbor]).norm();
            const float next_dist = curr_dist + edge_length;
            if (next_dist < distance[neighbor]) {
                distance[neighbor] = next_dist;
                pq.push({next_dist, neighbor});
            }
        }
    }

    return distance;
}

std::vector<int> assign_faces_to_sources(const std::vector<Vector3f>& centroids,
                                         const std::vector<std::vector<int>>& adjacency,
                                         const std::vector<int>& source_faces)
{
    const int num_faces = static_cast<int>(centroids.size());
    const float inf = std::numeric_limits<float>::infinity();

    std::vector<float> distance(num_faces, inf);
    std::vector<int> owner(num_faces, -1);

    struct QueueEntry {
        float distance = 0.0f;
        int face = -1;
        int source = -1;
    };

    struct CompareEntry {
        bool operator()(const QueueEntry& a, const QueueEntry& b) const
        {
            if (a.distance != b.distance) {
                return a.distance > b.distance;
            }
            return a.source > b.source;
        }
    };

    std::priority_queue<QueueEntry, std::vector<QueueEntry>, CompareEntry> pq;
    for (int source = 0; source < static_cast<int>(source_faces.size()); source++) {
        const int face = source_faces[source];
        distance[face] = 0.0f;
        owner[face] = source;
        pq.push({0.0f, face, source});
    }

    while (!pq.empty()) {
        const QueueEntry curr = pq.top();
        pq.pop();

        if (curr.distance > distance[curr.face]) {
            continue;
        }
        if (curr.distance == distance[curr.face] && curr.source != owner[curr.face]) {
            continue;
        }

        for (int neighbor : adjacency[curr.face]) {
            const float edge_length = (centroids[curr.face] - centroids[neighbor]).norm();
            const float next_dist = curr.distance + edge_length;

            if (next_dist < distance[neighbor] ||
                (std::abs(next_dist - distance[neighbor]) < 1e-6f && curr.source < owner[neighbor])) {
                distance[neighbor] = next_dist;
                owner[neighbor] = curr.source;
                pq.push({next_dist, neighbor, curr.source});
            }
        }
    }

    return owner;
}

float quantile_value(const std::vector<float>& sorted_values, float q)
{
    if (sorted_values.empty()) {
        return 0.0f;
    }

    const float clamped_q = std::clamp(q, 0.0f, 1.0f);
    const float idx = clamped_q * static_cast<float>(sorted_values.size() - 1);
    const int lo = static_cast<int>(std::floor(idx));
    const int hi = static_cast<int>(std::ceil(idx));
    const float t = idx - static_cast<float>(lo);

    return (1.0f - t) * sorted_values[lo] + t * sorted_values[hi];
}

std::vector<FaceComponent> extract_components(const std::vector<char>& active,
                                              const std::vector<std::vector<int>>& adjacency,
                                              const std::vector<Vector3f>& centroids,
                                              const std::vector<float>& distances,
                                              int min_component_size)
{
    std::vector<FaceComponent> components;
    std::vector<char> visited(active.size(), false);

    for (int start = 0; start < static_cast<int>(active.size()); start++) {
        if (!active[start] || visited[start]) {
            continue;
        }

        std::queue<int> q;
        q.push(start);
        visited[start] = true;

        FaceComponent component;

        while (!q.empty()) {
            const int face = q.front();
            q.pop();

            component.faces.push_back(face);
            component.centroid += centroids[face];
            component.max_distance = std::max(component.max_distance, distances[face]);
            component.mean_distance += distances[face];

            for (int neighbor : adjacency[face]) {
                if (active[neighbor] && !visited[neighbor]) {
                    visited[neighbor] = true;
                    q.push(neighbor);
                }
            }
        }

        if (static_cast<int>(component.faces.size()) < min_component_size) {
            continue;
        }

        const float inv_size = 1.0f / static_cast<float>(component.faces.size());
        component.centroid *= inv_size;
        component.mean_distance *= inv_size;
        component.score = component.max_distance * (1.0f + 0.05f * std::sqrt(static_cast<float>(component.faces.size())));
        components.push_back(component);
    }

    return components;
}

std::vector<int> bbox_fallback_partition(const std::vector<Vector3f>& vertices,
                                         const std::vector<Vector3i>& faces)
{
    Vector3f bbox_min = vertices[0];
    Vector3f bbox_max = vertices[0];
    for (const auto& v : vertices) {
        bbox_min = bbox_min.cwiseMin(v);
        bbox_max = bbox_max.cwiseMax(v);
    }

    const Vector3f bbox_size = (bbox_max - bbox_min).cwiseMax(Vector3f::Constant(1e-6f));
    const int front_axis = (bbox_size[0] > bbox_size[2]) ? 0 : 2;
    const int side_axis = (front_axis == 0) ? 2 : 0;

    std::vector<int> region(faces.size(), REGION_SHELL);
    for (int f = 0; f < static_cast<int>(faces.size()); f++) {
        const Vector3f c = (vertices[faces[f][0]] +
                            vertices[faces[f][1]] +
                            vertices[faces[f][2]]) / 3.0f;
        const Vector3f t = (c - bbox_min).cwiseQuotient(bbox_size);

        const float front = 1.0f - t[front_axis];
        const float side = t[side_axis];

        if (front > 0.85f && side > 0.35f && side < 0.65f) {
            region[f] = REGION_HEAD;
        } else if (front > 0.30f && side < 0.20f) {
            region[f] = REGION_LEFT_FRONT_FLIPPER;
        } else if (front > 0.30f && side > 0.80f) {
            region[f] = REGION_RIGHT_FRONT_FLIPPER;
        } else if (front < 0.30f && side < 0.20f) {
            region[f] = REGION_LEFT_REAR_FLIPPER;
        } else if (front < 0.30f && side > 0.80f) {
            region[f] = REGION_RIGHT_REAR_FLIPPER;
        }
    }

    return region;
}

Vector3f compute_axis_from_shell(const std::vector<Vector3f>& samples,
                                 bool major_axis)
{
    if (samples.size() < 3) {
        return Vector3f::Zero();
    }

    Vector2f mean = Vector2f::Zero();
    for (const auto& s : samples) {
        mean += Vector2f(s.x(), s.z());
    }
    mean /= static_cast<float>(samples.size());

    Matrix2f cov = Matrix2f::Zero();
    for (const auto& s : samples) {
        const Vector2f d = Vector2f(s.x(), s.z()) - mean;
        cov += d * d.transpose();
    }

    SelfAdjointEigenSolver<Matrix2f> solver(cov);
    if (solver.info() != Success) {
        return Vector3f::Zero();
    }

    const int col = major_axis ? 1 : 0;
    Vector2f axis2 = solver.eigenvectors().col(col);
    if (axis2.norm() < 1e-6f) {
        return Vector3f::Zero();
    }

    axis2.normalize();
    return Vector3f(axis2.x(), 0.0f, axis2.y());
}

void grow_head_region(std::vector<int>& region,
                      const std::vector<std::vector<int>>& adjacency,
                      const std::vector<Vector3f>& face_centroids,
                      const Vector3f& shell_center,
                      const Vector3f& front_axis,
                      const Vector3f& side_axis)
{
    std::vector<float> front_values;
    std::vector<float> side_values;
    front_values.reserve(face_centroids.size());
    side_values.reserve(face_centroids.size());

    float front_min = std::numeric_limits<float>::max();
    float front_max = -std::numeric_limits<float>::max();
    float head_front_sum = 0.0f;
    int head_count = 0;

    for (int f = 0; f < static_cast<int>(face_centroids.size()); f++) {
        const Vector3f offset = face_centroids[f] - shell_center;
        const float front = offset.dot(front_axis);
        const float side = std::abs(offset.dot(side_axis));

        front_values.push_back(front);
        side_values.push_back(side);
        front_min = std::min(front_min, front);
        front_max = std::max(front_max, front);

        if (region[f] == REGION_HEAD) {
            head_front_sum += front;
            head_count++;
        }
    }

    if (head_count == 0) {
        return;
    }

    std::sort(front_values.begin(), front_values.end());
    std::sort(side_values.begin(), side_values.end());

    const float front_range = std::max(1e-6f, front_max - front_min);
    const float head_front_mean = head_front_sum / static_cast<float>(head_count);
    const float front_threshold = std::max(
        quantile_value(front_values, 0.62f),
        head_front_mean - 0.28f * front_range);
    const float side_threshold = std::max(
        quantile_value(side_values, 0.28f),
        0.28f * (side_values.empty() ? 0.0f : side_values.back()));

    std::queue<int> q;
    std::vector<char> visited(region.size(), false);
    for (int f = 0; f < static_cast<int>(region.size()); f++) {
        if (region[f] == REGION_HEAD) {
            q.push(f);
            visited[f] = true;
        }
    }

    while (!q.empty()) {
        const int face = q.front();
        q.pop();

        for (int neighbor : adjacency[face]) {
            if (visited[neighbor]) {
                continue;
            }
            visited[neighbor] = true;

            if (region[neighbor] != REGION_SHELL) {
                continue;
            }

            const Vector3f offset = face_centroids[neighbor] - shell_center;
            const float front = offset.dot(front_axis);
            const float side = std::abs(offset.dot(side_axis));

            if (front >= front_threshold && side <= side_threshold) {
                region[neighbor] = REGION_HEAD;
                q.push(neighbor);
            }
        }
    }
}

void reinforce_head_band(std::vector<int>& region,
                         const std::vector<Vector3f>& vertices,
                         const std::vector<Vector3i>& faces)
{
    if (region.size() != faces.size()) {
        return;
    }

    const std::vector<Vector3f> face_centroids = compute_face_centroids(vertices, faces);
    const std::vector<std::vector<int>> adjacency = build_face_adjacency(faces);

    std::vector<Vector3f> shell_samples;
    Vector3f shell_center = Vector3f::Zero();
    Vector3f head_center = Vector3f::Zero();
    int head_count = 0;

    for (int f = 0; f < static_cast<int>(faces.size()); f++) {
        if (region[f] == REGION_SHELL || region[f] == REGION_HEAD) {
            shell_samples.push_back(face_centroids[f]);
            shell_center += face_centroids[f];
        }
        if (region[f] == REGION_HEAD) {
            head_center += face_centroids[f];
            head_count++;
        }
    }

    if (shell_samples.empty() || head_count == 0) {
        return;
    }

    shell_center /= static_cast<float>(shell_samples.size());
    head_center /= static_cast<float>(head_count);

    Vector3f front_axis = compute_axis_from_shell(shell_samples, true);
    if (front_axis.norm() < 1e-6f) {
        front_axis = Vector3f::UnitX();
    } else {
        front_axis.normalize();
    }

    Vector3f side_axis = compute_axis_from_shell(shell_samples, false);
    if (side_axis.norm() < 1e-6f) {
        side_axis = Vector3f(-front_axis.z(), 0.0f, front_axis.x());
    } else {
        side_axis.normalize();
    }

    if ((head_center - shell_center).dot(front_axis) < 0.0f) {
        front_axis = -front_axis;
    }

    std::vector<float> front_values;
    std::vector<float> side_values;
    front_values.reserve(faces.size());
    side_values.reserve(faces.size());
    for (int f = 0; f < static_cast<int>(faces.size()); f++) {
        const Vector3f offset = face_centroids[f] - shell_center;
        front_values.push_back(offset.dot(front_axis));
        side_values.push_back(std::abs(offset.dot(side_axis)));
    }

    std::vector<float> sorted_front = front_values;
    std::vector<float> sorted_side = side_values;
    std::sort(sorted_front.begin(), sorted_front.end());
    std::sort(sorted_side.begin(), sorted_side.end());

    const float front_threshold = quantile_value(sorted_front, 0.78f);
    const float side_threshold = quantile_value(sorted_side, 0.20f);

    for (int f = 0; f < static_cast<int>(faces.size()); f++) {
        if (region[f] == REGION_SHELL &&
            front_values[f] >= front_threshold &&
            side_values[f] <= side_threshold) {
            region[f] = REGION_HEAD;
        }
    }
}

std::vector<int> partition_single_pose(const std::vector<Vector3f>& vertices,
                                       const std::vector<Vector3i>& faces)
{
    if (vertices.empty() || faces.empty()) {
        return {};
    }

    const std::vector<Vector3f> face_centroids = compute_face_centroids(vertices, faces);
    const std::vector<std::vector<int>> adjacency = build_face_adjacency(faces);

    const Vector3f mesh_centroid = compute_mesh_centroid(vertices);
    int center_face = 0;
    float best_center_dist = std::numeric_limits<float>::max();
    for (int f = 0; f < static_cast<int>(face_centroids.size()); f++) {
        const float d = (face_centroids[f] - mesh_centroid).squaredNorm();
        if (d < best_center_dist) {
            best_center_dist = d;
            center_face = f;
        }
    }

    const std::vector<float> face_distances = compute_face_distances(face_centroids, adjacency, center_face);
    std::vector<float> sorted_distances = face_distances;
    std::sort(sorted_distances.begin(), sorted_distances.end());
    const float shell_expand_threshold = quantile_value(sorted_distances, 0.45f);

    const int min_component_size = std::max(10, static_cast<int>(faces.size() / 400));

    bool found_candidate = false;
    std::vector<FaceComponent> tip_components;

    for (float q = 0.95f; q >= 0.72f; q -= 0.02f) {
        const float threshold = quantile_value(sorted_distances, q);
        std::vector<char> active(faces.size(), false);
        for (int f = 0; f < static_cast<int>(faces.size()); f++) {
            active[f] = face_distances[f] >= threshold;
        }

        std::vector<FaceComponent> components = extract_components(
            active, adjacency, face_centroids, face_distances, min_component_size);

        if (static_cast<int>(components.size()) >= 5) {
            found_candidate = true;
            tip_components = std::move(components);
            break;
        }
    }

    if (!found_candidate) {
        return bbox_fallback_partition(vertices, faces);
    }

    std::sort(tip_components.begin(), tip_components.end(),
              [](const FaceComponent& a, const FaceComponent& b) {
                  return a.score > b.score;
              });

    if (tip_components.size() > 5) {
        tip_components.resize(5);
    }

    std::vector<int> source_faces;
    source_faces.reserve(6);
    source_faces.push_back(center_face);

    for (const auto& component : tip_components) {
        int seed_face = component.faces[0];
        float max_distance = -std::numeric_limits<float>::max();
        for (int face : component.faces) {
            if (face_distances[face] > max_distance) {
                max_distance = face_distances[face];
                seed_face = face;
            }
        }
        source_faces.push_back(seed_face);
    }

    if (static_cast<int>(source_faces.size()) != 6) {
        return bbox_fallback_partition(vertices, faces);
    }

    const std::vector<int> source_owner = assign_faces_to_sources(face_centroids, adjacency, source_faces);
    std::vector<int> region(faces.size(), REGION_SHELL);

    std::vector<Vector3f> shell_samples;
    shell_samples.reserve(faces.size());
    Vector3f shell_center = Vector3f::Zero();
    for (int f = 0; f < static_cast<int>(faces.size()); f++) {
        if (source_owner[f] == 0) {
            shell_samples.push_back(face_centroids[f]);
            shell_center += face_centroids[f];
        }
    }

    if (!shell_samples.empty()) {
        shell_center /= static_cast<float>(shell_samples.size());
    } else {
        shell_center = mesh_centroid;
        shell_samples = face_centroids;
    }

    Vector3f front_axis = compute_axis_from_shell(shell_samples, true);
    if (front_axis.norm() < 1e-6f) {
        front_axis = Vector3f::UnitX();
    } else {
        front_axis.normalize();
    }

    Vector3f side_axis = compute_axis_from_shell(shell_samples, false);
    if (side_axis.norm() < 1e-6f) {
        side_axis = Vector3f(-front_axis.z(), 0.0f, front_axis.x());
    } else {
        side_axis.normalize();
    }

    std::array<FaceComponent, 5> appendage_regions;
    for (int source = 1; source <= 5; source++) {
        for (int f = 0; f < static_cast<int>(faces.size()); f++) {
            if (source_owner[f] == source) {
                appendage_regions[source - 1].faces.push_back(f);
                appendage_regions[source - 1].centroid += face_centroids[f];
            }
        }

        if (appendage_regions[source - 1].faces.empty()) {
            return bbox_fallback_partition(vertices, faces);
        }

        appendage_regions[source - 1].centroid /=
            static_cast<float>(appendage_regions[source - 1].faces.size());
    }

    int head_component = 0;
    float head_score = std::numeric_limits<float>::max();
    for (int c = 0; c < static_cast<int>(appendage_regions.size()); c++) {
        const Vector3f offset = appendage_regions[c].centroid - shell_center;
        const float forward = std::abs(offset.dot(front_axis));
        const float side = std::abs(offset.dot(side_axis));
        const float score = side / (forward + 1e-4f);
        if (score < head_score) {
            head_score = score;
            head_component = c;
        }
    }

    if ((appendage_regions[head_component].centroid - shell_center).dot(front_axis) < 0.0f) {
        front_axis = -front_axis;
    }

    const Vector3f up_axis = Vector3f::UnitY();
    Vector3f right_axis = up_axis.cross(front_axis);
    if (right_axis.squaredNorm() > 1e-8f) {
        right_axis.normalize();
    }
    if (right_axis.squaredNorm() > 1e-8f && side_axis.dot(right_axis) < 0.0f) {
        side_axis = -side_axis;
    }

    struct RankedComponent {
        int index = -1;
        float front = 0.0f;
        float side = 0.0f;
    };

    std::array<int, 5> component_to_region = {
        REGION_SHELL,
        REGION_SHELL,
        REGION_SHELL,
        REGION_SHELL,
        REGION_SHELL
    };
    component_to_region[head_component] = REGION_HEAD;

    std::vector<RankedComponent> flipper_components;
    flipper_components.reserve(4);
    for (int c = 0; c < static_cast<int>(appendage_regions.size()); c++) {
        if (c == head_component) {
            continue;
        }

        const Vector3f offset = appendage_regions[c].centroid - shell_center;
        flipper_components.push_back({
            c,
            offset.dot(front_axis),
            offset.dot(side_axis)
        });
    }

    if (flipper_components.size() != 4) {
        return bbox_fallback_partition(vertices, faces);
    }

    std::sort(flipper_components.begin(), flipper_components.end(),
              [](const RankedComponent& a, const RankedComponent& b) {
                  return a.front > b.front;
              });

    std::array<RankedComponent, 2> front_pair = {
        flipper_components[0],
        flipper_components[1]
    };
    std::array<RankedComponent, 2> rear_pair = {
        flipper_components[2],
        flipper_components[3]
    };

    std::sort(front_pair.begin(), front_pair.end(),
              [](const RankedComponent& a, const RankedComponent& b) {
                  return a.side < b.side;
              });
    std::sort(rear_pair.begin(), rear_pair.end(),
              [](const RankedComponent& a, const RankedComponent& b) {
                  return a.side < b.side;
              });

    component_to_region[front_pair[0].index] = REGION_LEFT_FRONT_FLIPPER;
    component_to_region[front_pair[1].index] = REGION_RIGHT_FRONT_FLIPPER;
    component_to_region[rear_pair[0].index] = REGION_LEFT_REAR_FLIPPER;
    component_to_region[rear_pair[1].index] = REGION_RIGHT_REAR_FLIPPER;

    for (int c = 0; c < static_cast<int>(appendage_regions.size()); c++) {
        const int mapped_region = component_to_region[c];
        for (int face : appendage_regions[c].faces) {
            region[face] = mapped_region;
        }
    }

    grow_head_region(region, adjacency, face_centroids, shell_center, front_axis, side_axis);

    for (int iter = 0; iter < 3; iter++) {
        std::vector<int> next_region = region;

        for (int f = 0; f < static_cast<int>(faces.size()); f++) {
            if (region[f] != REGION_SHELL || face_distances[f] < shell_expand_threshold) {
                continue;
            }

            int neighbor_counts[NUM_REGIONS] = {0};
            for (int neighbor : adjacency[f]) {
                const int r = region[neighbor];
                if (r >= 0 && r < NUM_REGIONS && r != REGION_SHELL) {
                    neighbor_counts[r]++;
                }
            }

            int best_region = REGION_SHELL;
            int best_count = 0;
            for (int r = 0; r < REGION_SHELL; r++) {
                if (neighbor_counts[r] > best_count) {
                    best_count = neighbor_counts[r];
                    best_region = r;
                }
            }

            if (best_count >= 2) {
                next_region[f] = best_region;
            }
        }

        region.swap(next_region);
    }

    return region;
}

float partition_score(const std::vector<int>& region)
{
    if (region.empty()) {
        return -std::numeric_limits<float>::infinity();
    }

    int counts[NUM_REGIONS] = {0};
    for (int r : region) {
        if (r >= 0 && r < NUM_REGIONS) {
            counts[r]++;
        }
    }

    if (counts[REGION_HEAD] == 0 ||
        counts[REGION_LEFT_FRONT_FLIPPER] == 0 ||
        counts[REGION_RIGHT_FRONT_FLIPPER] == 0 ||
        counts[REGION_LEFT_REAR_FLIPPER] == 0 ||
        counts[REGION_RIGHT_REAR_FLIPPER] == 0) {
        return -std::numeric_limits<float>::infinity();
    }

    const int flipper_counts[4] = {
        counts[REGION_LEFT_FRONT_FLIPPER],
        counts[REGION_RIGHT_FRONT_FLIPPER],
        counts[REGION_LEFT_REAR_FLIPPER],
        counts[REGION_RIGHT_REAR_FLIPPER]
    };

    const int min_flipper = *std::min_element(std::begin(flipper_counts), std::end(flipper_counts));
    const int max_flipper = *std::max_element(std::begin(flipper_counts), std::end(flipper_counts));

    return static_cast<float>(min_flipper) -
           0.35f * static_cast<float>(max_flipper - min_flipper) +
           0.20f * static_cast<float>(std::min(counts[REGION_HEAD], 600));
}

} // namespace

std::string partition_labels_path()
{
    return "turtle_poses/frame_20_partition_labels.txt";
}

std::vector<int> load_partition_labels(const std::string& path, int expected_faces)
{
    std::ifstream in(path);
    if (!in.is_open()) {
        return {};
    }

    std::vector<int> region;
    region.reserve(expected_faces);

    int label = REGION_SHELL;
    while (in >> label) {
        if (label < 0 || label >= NUM_REGIONS) {
            return {};
        }
        region.push_back(label);
    }

    if (static_cast<int>(region.size()) != expected_faces) {
        return {};
    }

    return region;
}

bool save_partition_labels(const std::string& path, const std::vector<int>& region)
{
    std::filesystem::path out_path(path);
    std::error_code ec;
    if (!out_path.parent_path().empty()) {
        std::filesystem::create_directories(out_path.parent_path(), ec);
    }

    std::ofstream out(path);
    if (!out.is_open()) {
        return false;
    }

    for (int label : region) {
        out << label << '\n';
    }

    return out.good();
}

std::vector<int> partition_mesh(const std::vector<Vector3f>& vertices,
                                const std::vector<Vector3i>& faces)
{
    constexpr int kReferenceFrame = 20;

    struct CachedPartition {
        int face_count = -1;
        std::vector<int> region;
    };

    static CachedPartition cache;

    if (cache.face_count == static_cast<int>(faces.size()) &&
        cache.region.size() == faces.size()) {
        return cache.region;
    }

    const std::vector<int> saved_region = load_partition_labels(partition_labels_path(),
                                                                static_cast<int>(faces.size()));
    if (!saved_region.empty()) {
        cache.face_count = static_cast<int>(faces.size());
        cache.region = saved_region;
        return cache.region;
    }

    std::vector<int> best_region;
    const std::string reference_path = "turtle_poses/frame_" + std::to_string(kReferenceFrame) + ".obj";
    TriMesh reference_mesh = load_obj(reference_path);

    if (!reference_mesh.vertices.empty() && reference_mesh.faces.size() == faces.size()) {
        best_region = partition_single_pose(reference_mesh.vertices, reference_mesh.faces);
        reinforce_head_band(best_region, reference_mesh.vertices, reference_mesh.faces);
    } else {
        best_region = partition_single_pose(vertices, faces);
        reinforce_head_band(best_region, vertices, faces);
    }

    cache.face_count = static_cast<int>(faces.size());
    cache.region = best_region;
    return cache.region;
}
