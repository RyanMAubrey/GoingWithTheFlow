// h = time step
// gamma_k = all vertex positions at time k
// gamma_k1 = all vertex positions at time k+1
// faces = each face has i,j,k
// edges = edge struct

struct Edge {
    int i;
    int j;
    float alpha; // [-pi, pi]
};

vec6 calc_fluid_momentum(std::vector<vec3>& gamma_k, std::vector<vec3>& gamma_k1,
                         std::vector<ivec3>& faces, std::vector<Edge>& edges,
                         float rho_f, float h) {
    // Calculate face values
    std::vector<float> face_areas(faces.size());
    std::vector<vec3> face_normals(faces.size());
    for (int f = 0; f < faces.size(); f++) {
        int i = faces[f].x;
        int j = faces[f].y;
        int k = faces[f].z;
        vec3 normal = cross(gamma_k[j] - gamma_k[i], gamma_k[k] - gamma_k[i]);

        face_areas[f] = 0.5f * length(normal); // area = 1/2 * ||N||
        face_normals[f] = normalize(normal);
    }

    // Calculate delta
    float total_area = 0.0f;
    for (int f = 0; f < faces.size(); f++) {
        total_area += face_areas[f];
    }
    float total_edge = 0.0f;
    for (int e = 0; e < edges.size(); e++) {
        // Given edge index, find vertices for length & bending angle
        Edge curr_e = edges[e];
        float alpha_ij = curr_e.alpha;
        float l_ij = length(gamma_k[curr_e].j] - gamma_k[curr_e.i]);
        total_edge += alpha_ij * l_ij;
    }
    float delta = total_area / total_edge;

    vec3 l_f(0.0f);
    vec3 p_f(0.0f);
    for (int f = 0; f < faces.size(); f++) {
        int i = faces[f].x;
        int j = faces[f].y;
        int k = faces[f].z;
        // Combine ijk into one value (centers)
        vec3 gamma_ijk = (gamma_k[i] + gamma_k[j] + gamma_k[k]) / 3.0f;
        vec3 gamma_prime_ijk = ((gamma_k1[i] - gamma_k[i]) + (gamma_k1[j] - gamma_k[j]) + (gamma_k1[k] - gamma_k[k])) / (3.0f * h);

        vec3 n_ijk = face_normals[f];
        float A_ijk = face_areas[f];
        float dot_prod = dot(gamma_prime_ijk, n_ijk);

        l_f += rho_f * delta * dot_prod * cross(gamma_ijk, n_ijk) * A_ijk;
        p_f += rho_f * delta * dot_prod * n_ijk * A_ijk;
    }

    return vec6(l_f, p_f);
}
