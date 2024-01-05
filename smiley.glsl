float _smooth(float a, float b, float t) {
    return clamp((t-a)/(b-a), 0., 1.);
}

float _highlight(float a, float b, float c, float d, float t) {
    return clamp(((t-a)/(b-a))*(d-c) + c, 0., 1.);
}

vec2 within(vec2 uv, vec4 rect) {
    return (uv-rect.xy) / (rect.zw-rect.xy);
}

vec4 mouth(vec2 uv) {
    uv -= 0.5;
    vec4 col = vec4(.5, .18, .05, 1.);

    uv.y *= 1.5;
    uv.y -= uv.x * uv.x * 1.3;
    float d = length(uv);

    float td = length(uv-vec2(0., .6));
    vec3 toothCol = vec3(1.)  * smoothstep(.6, .35, d);
    col.rgb = mix(col.rgb, toothCol, smoothstep(.4, .37, td));

    td = length(uv+vec2(0., .5));
    col.rgb = mix(col.rgb, vec3(1., .5, .5), smoothstep(.5, .2, td));
    col.a = smoothstep(.5, .48, d);
    return col;
}

vec4 eye(vec2 uv) {
    uv -= 0.5;
    float d = length(uv);

    vec4 irisColor = vec4(.3, .5, 1., 1.);
    vec4 col = vec4(1.);
    col = mix(col, irisColor, smoothstep(.1, .7, d) * .5);
    
    col.rgb *= 1. - smoothstep(.45, .5, d) * .5 * clamp(-uv.y-uv.x, 0., 1.);
    col.rgb = mix(col.rgb, vec3(0.), smoothstep(.3, .28, d));
    irisColor.rgb *= 1. + smoothstep(.3, .05, d);
    col.rgb = mix(col.rgb, irisColor.rgb, smoothstep(.28, .25, d));
    col.rgb = mix(col.rgb, vec3(0.), smoothstep(.16, .14, d));

    float highlight = smoothstep(.1, .09, length(uv-vec2(-.15, .15)));
    highlight += smoothstep(.07, .05, length(uv+vec2(-.08, .08)));
    col.rgb =  mix(col.rgb, vec3(1.), highlight);

    col.a = smoothstep(.5, .48, d);
    return col;
}

vec4 head(vec2 uv) {
    vec4 col = vec4(.9, .65, .1, 1.);
    float d = length(uv);
    col.a = smoothstep(.5, .49, d);

    float edgeshade = _smooth(.35, .5, d);
    edgeshade *= edgeshade;
    col.rgb *= 1. - edgeshade * .5;
    col.rgb = mix(col.rgb, vec3(.6, .3, .1), smoothstep(.47, .48, d));

    float highlight = smoothstep(.41, .405, d);
    highlight *= _highlight(.41, .0, .75, .0, uv.y);
    highlight *= smoothstep(.18, .19, length(uv-vec2(.21, .08)));
    col.rgb = mix(col.rgb, vec3(1.), highlight);

    d = length(uv - vec2(.25, -.2));
    float cheek = smoothstep(.2, .01, d) * .4;
    cheek *= smoothstep(.18, .17, d);
    col.rgb = mix(col.rgb, vec3(1., .1, .1), cheek);
    return col;
}

vec4 smiley(vec2 uv) {
    vec4 col = vec4(0.);
    uv.x = abs(uv.x);
    
    vec4 headColor = head(uv);
    col = mix(col, headColor, headColor.a);

    vec4 eyeColor = eye(within(uv, vec4(.03, -.1, .37, .25)));
    col = mix(col, eyeColor, eyeColor.a);

    vec4 mouthColor = mouth(within(uv, vec4(-.3, -.4, .3, -.1)));
    col = mix(col, mouthColor, mouthColor.a);
    
    return col;
}

void mainImage(out vec4 fragColor, in vec2 fragCoord) {
    vec2 uv = fragCoord / iResolution.xy - 0.5;
    float aspect = iResolution.x / iResolution.y;
    if (aspect < 1.) uv.y /= aspect; else uv.x *= aspect;
    
    fragColor = smiley(uv);
}