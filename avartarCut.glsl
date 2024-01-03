#iChannel0 "file://Textures/avatar.jpg"

void mainImage(out vec4 fragColor, in vec2 fragCoord) {
    vec2 uv = fragCoord.xy / iResolution.xy - 0.5;
    if (iResolution.x < iResolution.y) {
        uv.y *= iResolution.y / iResolution.x;
    } else {
        uv.x *= iResolution.x / iResolution.y;
    }
    fragColor = texture2D(iChannel0, uv + 0.5);
    fragColor *= smoothstep(0.5, 0.45, length(uv));
}