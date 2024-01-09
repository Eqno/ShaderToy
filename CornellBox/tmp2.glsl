float DistributionGGX(const float dnh, const float r) {
    float a = r * r, _dnh = max(0., dnh), dnh_square = _dnh * _dnh;
    return a * a / max(PI * pow(dnh_square * (a * a - 1.) + 1., 2.), EPS);
}

vec3 FresnelSchlick(const float dvh, const vec3 f) {
    return f + (1. - f) * pow(saturate(1. - dvh), 5.);
}

float GeometrySchlickGGX(float dnx, float k) {
    return dnx / max(dnx * (1. - k) + k, EPS);
}

float GeometrySmith(float dnv, float dnl, float r) {
    return GeometrySchlickGGX(saturate(dnv), r / 2.)
        * GeometrySchlickGGX(saturate(dnl), r / 2.);
}

vec3 BRDF(const vec3 n, const vec3 l, const vec3 v, Material mat)
{
    vec3 h = normalize(v + l);
    float dvh = dot(v, h);
    float dnh = dot(n, h);
    float dnl = dot(n, l);
    float dnv = dot(n, v);
    vec3 fresnel = mix(vec3(.4), mat.color, mat.metallic);

    float D = DistributionGGX(dnh, mat.roughness);
    vec3 F = FresnelSchlick(max(dvh, 0.), fresnel);
    float G = GeometrySmith(dnv, dnl, mat.roughness);

    vec3 specular = D * G * F / max(4. * saturate(dnv) * saturate(dnl), EPS);
    vec3 kd = mix(vec3(1.) - F, vec3(0.), mat.metallic);
    return kd * mat.color / PI + specular;
}
