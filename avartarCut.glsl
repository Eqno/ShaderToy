#iChannel0 "file://Textures/avatar.jpg"

void mainImage(out vec4 fragColor, in vec2 fragCoord) {
    vec2 uv = fragCoord / iResolution.xy - 0.5;
    float aspect = iResolution.x / iResolution.y;
    if (aspect < 1.) uv.y /= aspect; else uv.x *= aspect;
    
    fragColor = texture2D(iChannel0, uv + 0.5);
    fragColor *= smoothstep(0.5, 0.45, length(uv));
}