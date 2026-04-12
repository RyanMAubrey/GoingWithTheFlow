// h = time step
// gamma_k = all vertex positions at time k
// gamma_k1 = all vertex positions at time k+1
// mass density = non-negative

vec6 calc_body_momentum(std::vector<vec3>& gamma_k, std::vector<vec3>& gamma_k1, std::vector<float>& mass_density, float h) {
    vec3 l_b(0.0f);
    vec3 p_b(0.0f);

    for (int i = 0; i < gamma_k.size(); i++) {
        float rho_b = mass_density[i];
        vec3 gamma_prime_i = (gamma_k1[i] - gamma_k[i]) / h;

        l_b += cross(gamma_k[i], gamma_prime_i) * rho_b;
        p_b += gamma_prime_i * rho_b;
    }

    return vec6(l_b, p_b);
}
