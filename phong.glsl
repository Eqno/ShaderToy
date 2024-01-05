vec3 background(vec3 rd) {
    float sky = max(0., dot(rd, vec3(0., 1., 0.)));
    vec3 col = vec3(.5, .8, 1.);
    return pow(sky, 1.) * col;
}

void mainImage(out vec4 fragColor, in vec2 fragCoord) {
    vec2 uv = fragCoord / iResolution.xy - .5;
    float aspect = iResolution.x / iResolution.y;
    if (aspect < 1.) uv.y /= aspect; else uv.x *= aspect;

    vec3 col = vec3(.6);
    vec3 camera = vec3(0., 0., 0);
    vec3 center = vec3(0., 0., 5.);
    float radius = .3;

    if (length(uv-center.xy) <= radius) {
        col = vec3(.8, .3, 1);
        vec3 pointlight = vec3(10.*sin(iTime*.5), 5., -1.);
        float z = center.z - sqrt(pow(radius, 2.) - pow(uv.y-center.y, 2.));

        vec3 curr = vec3(uv, z);
        vec3 normal = normalize(curr - center);
        vec3 light = normalize(pointlight - curr);
        vec3 view = normalize(camera - curr);

        vec3 h =  normalize(light + view);
        vec3 ambient = .1 * col;
        vec3 diffuse = vec3(max(dot(normal, light), 0.));
        vec3 specular = vec3(max(.2 * pow(dot(normal, h), 200.), 0.));
        col *= ambient + diffuse + specular;
    } else {
        col = background(vec3(uv, 0.));
    }
    fragColor = vec4(col, 1.);
}