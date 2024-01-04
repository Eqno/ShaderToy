#define ANTIALIAS_ALWAYS true
#define ANTIALIAS_SAMPLES 4
#define ANIMATE_CAMERA true

#define PI 3.14159265359
#define EXPOSURE 34.
#define GAMMA 2.1
#define SOFT_SHADOWS_FACTOR 4.
#define MAX_RAYMARCH_ITER 128
#define MAX_RAYMARCH_ITER_SHADOWS 16
#define MIN_RAYMARCH_DELTA 0.0015
#define GRADIENT_DELTA 0.0002
#define OBJ_FLOOR		1.
#define OBJ_CEILING		2.
#define OBJ_BACKWALL	3.
#define OBJ_LEFTWALL	4.
#define OBJ_RIGHTWALL	5.
#define OBJ_LIGHT		6.
#define OBJ_SHORTBLOCK	7.
#define OBJ_TALLBLOCK	8.

const vec3 lightColor = vec3(16.86, 8.76 +2., 3.2 + .5);
const vec3 lightDiffuseColor = vec3(.78);
const vec3 leftWallColor = vec3(.611, .0555, .062);
const vec3 rightWallColor = vec3(.117, .4125, .115);
const vec3 whiteWallColor = vec3(.7295, .7355, .729);
const vec3 cameraTarget = vec3(556, 548.8, 559.2) * .5;

float sd_box(vec3 p, vec3 b)  {
    vec3 d = abs(p) - b;
    return min(max(d.x, max(d.y, d.z)), 0.) + length(max(d, 0.));
}

vec3 rotateX(vec3 p, float a) {
    float c = cos(a), s = sin(a);
    return vec3(p.x, c*p.y-s*p.z, s*p.y+c*p.z);
}

vec3 rotateY(vec3 p, float a) {
    float c = cos(a), s = sin(a);
    return vec3(c*p.x+s*p.z, p.y, -s*p.x+c*p.z);
}

vec3 rotateZ(vec3 p, float a) {
    float c = cos(a), s = sin(a);
    return vec3(c*p.x-s*p.y, s*p.x+c*p.y, p.z);
}

vec2 map_blocks(vec3 p, vec3 ray_d) {
    vec2 res = vec2(OBJ_SHORTBLOCK, sd_box(rotateY(p + vec3(-186, -82.5, -169.5), 0.29718), vec3(83.66749, 83.522452, 82.5)));
	vec2 obj1 = vec2(OBJ_TALLBLOCK, sd_box(rotateY(p + vec3(-368.5, -165, -351.5), -0.30072115), vec3(87.02012, 165, 83.6675)));
	if (obj1.y < res.y) res = obj1;
	return res;
}

vec2 map(vec3 p, vec3 ray_d) {
    vec2 res = vec2(OBJ_FLOOR, p.y);
	vec2 obj1 = vec2(OBJ_CEILING, 548.8 - p.y);
	if (obj1.y < res.y) res = obj1;
	vec2 obj2 = vec2(OBJ_BACKWALL, 559.2 - p.z);
	if (obj2.y < res.y) res = obj2;
	vec2 obj3 = vec2(OBJ_LEFTWALL, 556. - p.x);
	if (obj3.y < res.y) res = obj3;
	vec2 obj4 = vec2(OBJ_RIGHTWALL, p.x);
	if (obj4.y < res.y) res = obj4;
	vec2 obj5 = vec2(OBJ_LIGHT, sd_box(p + vec3(-278, -548.8, -292), vec3(65, 0.05, 65)));
	if (obj5.y < res.y) res = obj5;
	vec2 obj6 = map_blocks(p, ray_d);
	if (obj6.y < res.y) res = obj6;
	return res;
}

vec2 map(vec3 p) {
    return map(p, vec3(0.));
}

vec3 gradientNormal(vec3 p) {
    return normalize(vec3(
        map(p+vec3(GRADIENT_DELTA, 0., 0.)).y-map(p-vec3(GRADIENT_DELTA, 0., 0.)).y, 
        map(p+vec3(0., GRADIENT_DELTA, 0.)).y-map(p-vec3(0., GRADIENT_DELTA, 0.)).y, 
        map(p+vec3(0., 0., GRADIENT_DELTA)).y-map(p-vec3(0., 0., GRADIENT_DELTA)).y));
}

float raymarch(vec3 ray_s, vec3 ray_d, out float d, out vec3 p, out int iter) {
    d = 0.;
    float min_step = .1;
    vec2 map_res;
    for (int i=1; i<MAX_RAYMARCH_ITER; i++) {
        p = ray_s + ray_d * d;
        map_res = map(p, ray_d);
        if (map_res.y < MIN_RAYMARCH_DELTA) {
            iter = i;
            return map_res.x;
        }
        d += max(map_res.y, min_step);
    }
    return -1.;
}

vec3 render(vec3 ray_s, vec3 ray_d) {
    float d;
    vec3 p;
    int iter;
    float object_id = raymarch(ray_s, ray_d, d, p, iter);

    vec3 col = vec3(0.);
    if (p.z >= 0.) {
        if (object_id == OBJ_FLOOR) col = whiteWallColor;
		else if (object_id == OBJ_CEILING) col = whiteWallColor;
		else if (object_id == OBJ_BACKWALL) col = whiteWallColor;
		else if (object_id == OBJ_LEFTWALL) col = leftWallColor;
		else if (object_id == OBJ_RIGHTWALL) col = rightWallColor;
		else if (object_id == OBJ_LIGHT) col = lightDiffuseColor;
		else if (object_id == OBJ_SHORTBLOCK) col = whiteWallColor;
		else if (object_id == OBJ_TALLBLOCK) col = whiteWallColor;

        if (object_id == OBJ_LIGHT) {
            col *= lightColor;
        } else {
            float light_size = 25.;
            vec3 light_pos = vec3(278, 548.8 -50., 292 - 50);
            if (object_id == OBJ_CEILING) light_pos.y -= 550.;

            light_pos.x = max(light_pos.x - light_size, min(light_pos.x + light_size, p.x));
            light_pos.y = max(light_pos.y - light_size, min(light_pos.y + light_size, p.y));
            vec3 n = gradientNormal(p);

            vec3 l =  normalize(light_pos - p);
            float light_dis = length(light_pos - p);
            float atten = ((1. / light_dis) * .5) + ((1. / (light_dis * light_dis)) * .5);

            
        }
    }

}

vec3 rotateCamera(vec3 ray_s, vec3 ray_d) {
    ray_d.x = -ray_d.x;
    vec3 target = normalize(cameraTarget - ray_s);

    float ang_y = atan(target.z, target.x);
    ray_d = rotateY(ray_d, PI / 2. - ang_y);
    
    float ang_x = atan(target.y, target.z);
    ray_d = rotateX(ray_d, -ang_x);

    if (ANIMATE_CAMERA) {
        float ang_z = smoothstep(0., 1., (iTime-5.)*.1) * sin(iTime*1.1+.77) * .05;
        ray_d = rotateZ(ray_d, ang_z);
    }
    return ray_d;
}

vec3 moveCamera(vec3 camera_start) {
    camera_start += vec3(278, 273, -800);
    if (ANIMATE_CAMERA) {
        vec3 ray_s = camera_start + vec3(cos(iTime*.8)*180., cos(iTime*.9)*180., (cos(iTime*.3)+1.)*390.);
        return mix(camera_start, ray_s, smoothstep(0., 1., (iTime-5.)*.1));
    } else return camera_start;
}

void mainImage(out vec4 fragColor, in vec2 fragCoord) {
    vec3 camera_start = vec3(0., 0., -1.4);
    vec3 col = vec3(0);
    if (ANTIALIAS_ALWAYS || iTime < 5.) {
        float d_ang =  2. * PI / float(ANTIALIAS_SAMPLES);
        float ang = d_ang * .333;
        float r = .4;
        for (int i=0; i<ANTIALIAS_SAMPLES; i++) {
            vec2 position = vec2((fragCoord.x+cos(ang)*r-iResolution.x*.5)/iResolution.y, (fragCoord.y+sin(ang)*r-iResolution.y*.5)/iResolution.y);
            vec3 ray_s = moveCamera(camera_start);
            vec3 ray_d = rotateCamera(ray_s, normalize(vec3(position, 0) - camera_start));
            col += render(ray_s, ray_d);
            ang += d_ang;
        }
        col /= float(ANTIALIAS_SAMPLES);
    }
    col *= EXPOSURE;
    col = pow(col, vec3(1. / GAMMA));
    fragColor = vec4(col, 1.);
}