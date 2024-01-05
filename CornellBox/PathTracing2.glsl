// 数学
#define PI 3.141592653589
#define INFINITY 9999999.
#define RANDOM_SEED 0

struct Random {
    int seed;
    int flat_index;
} random;

void updateRandom() {
    random.seed = 0;
    random.flat_index = int(dot(gl_FragCoord.xy, vec2(1., 4096.)));
}

void EncryptTea(inout uvec2 arg) {
    uvec4 key = uvec4(0xa341316c, 0xc8013ea4, 0xad90777d, 0x7e95761e);
    uint v0 = arg[0], v1 = arg[1], sum = 0u, delta = 0x9e3779b9u;
    for (int i=0; i<32; i++) {
        sum += delta;
        v0 += ((v1<<4)+key[0]) ^ (v1+sum) ^ ((v1>>5)+key[1]);
        v1 += ((v0<<4)+key[2]) ^ (v0+sum) ^ ((v0>>5)+key[3]);
    }
    arg[0] = v0;
    arg[1] = v1;
}

vec2 GetRandom() {
    uvec2 arg = uvec2(random.flat_index, random.seed++);
    EncryptTea(arg);
    return fract(vec2(arg) / vec2(0xffffffffu));
}

// 光照
#define LIGHT_SIZE .5;
#define LIGHT_NORMAL vec3(0., -1., 0.)
#define LIGHT_POSITION vec3(.5 * sin(iTime), .9, .5 * cos(iTime))

struct Light {
    float size;
    float area;
    vec3 normal;
    vec4 albedo;
    vec3 position;
} light;

void updateLight() {
#ifdef LIGHT_SIZE
    light.size = LIGHT_SIZE;
#else
    light.size = .5;
#endif
#ifdef LIGHT_NORMAL
    light.normal = LIGHT_NORMAL;
#else
    light.normal = vec3(0., -1., 0.);
#endif
#ifdef LIGHT_POSITION
    light.position = LIGHT_POSITION;
#else
    light.position = vec3(0., .9, 0.);
#endif
    light.area = light.size * light.size;
    light.albedo = vec4(1., 1., 1., 2. / light.area);
}

// 追踪
struct Ray {
    vec3 origin, direction;
};

struct AABB {
    vec3 min_point, max_point;
};

bool intersectAABB(Ray ray, AABB aabb, inout float min_t, inout float max_t) {

    // 先求逆再乘，两次除法变一次除法 + 两次乘法，性能更高
    vec3 inv_dir = 1. / ray.direction;

    // 解 AABB 两个顶点到射线出发点，以射线的方向为基的 t
    vec3 t_1 = (aabb.min_point - ray.origin) * inv_dir;
	vec3 t_2 = (aabb.max_point - ray.origin) * inv_dir;

    // 分别对 t1、t2 每个分量取 min、max 重新组成 vec3
    vec3 t_min = min(t_1, t_2);
    vec3 t_max = max(t_1, t_2);

    // 对上一步结果分量再取 min、max，检查是否有数值交叉
    min_t = max(max(min_t, t_min.x), max(t_min.y, t_min.z));
    max_t = min(min(max_t, t_max.x), min(t_max.y, t_max.z));

    // 根据 slab 算法，若没有分量交叉，则射线与 AABB 相交
    return min_t < max_t;
}

float intersectPlane(Ray ray, vec3 point, vec3 normal) {

    // dot(origin + t * direction, normal) = dot(point - origin, normal) 
    float denom = dot(ray.direction, normal);
    float t = dot(point - ray.origin, normal) / denom;
    return t > 0. ? t : INFINITY;
}

// 渲染
#define CAMERA_POSITION vec3(0., 0., 3.125)
#define NUM_SAMPLES 128
#define NUM_BOUNCES 2

vec3 sampleColor(Ray ray) {
    vec3 contribution = vec3(0.);
    vec3 tp = vec3(1.);

    vec3 position, normal;
    vec4 albedo;
    float t = intersect(ray, position, normal, albedo);

    for (int i=0; i<NUM_BOUNCES; i++) {

        vec3 pos_ls = sampleLight(GetRandom());
        vec3 l_nee = pos_ls - position;
    }
}

// 此程序（Pixel Shader）的入口函数
void mainImage(out vec4 fragColor, in vec2 fragCoord) {

    // 为当前片元着色更新灯光和随机数据
    updateLight();
    updateRandom();

    // 将 uv 移动至屏幕中间并处理长宽比
    vec2 uv = fragCoord / iResolution.xy - .5;
    float aspect = iResolution.x / iResolution.y;
    if (aspect < 1.) uv.y /= aspect; else uv.x *= aspect;

    // 由此片元发出若干条射线进行路径追踪
    vec3 color = vec3(0.);
    for (int i=0; i<NUM_SAMPLES; i++) {

        // 从摄像机出发
        Ray ray;
        ray.origin = CAMERA_POSITION;
        
        // 随机方向投射
        vec2 random_result = GetRandom();
        ray.direction = normalize(vec3(
            uv+random_result.x*dFdx(uv) + random_result.y*dFdy(uv),
            -1.
        ));
        
        // 将每次追踪采样得到的颜色值累加
        color += sampleColor(ray);
    }

    // 算出平均颜色并做 Gamma 运算得到片元结果
    fragColor = vec4(pow(color/float(NUM_SAMPLES), vec3(1./2.2)), 1.);
}