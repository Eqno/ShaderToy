// 数学
#define PI 3.141592653589                       // 常量 π
#define INFINITY 9999999.                       // 无限大
#define RANDOM_SEED 0                           // 初始随机种子

struct Random {
    int seed;                                   // 随机种子
    int flat_index;                             // 随机索引
} random;

void updateRandom() {

    // 随机参数 rand seed 和 index
    random.seed = RANDOM_SEED;
    random.flat_index = int(dot(gl_FragCoord.xy, vec2(1., 4096.)));
}

void EncryptTea(inout uvec2 arg) {

    // 对输入的 arg 进行 TEA 加密
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

    // fract 取小数部分获得随机数
    uvec2 arg = uvec2(random.flat_index, random.seed++);
    EncryptTea(arg);
    return fract(vec2(arg) / vec2(0xffffffffu));
}

// 光照
struct Light {
    float size;                                 // 光源大小
    float area;                                 // 光源区域
    vec3 normal;                                // 光源朝向
    vec4 albedo;                                // 光反射率
    vec3 position;                              // 光源位置
} light;

void updateLight() {
    light.size = .5;
    light.area = light.size * light.size;
    light.normal = vec3(0., -1., 0.);
    light.albedo = vec4(1., 1., 1., 2. / light.area);
    light.position = vec3(.5 * sin(iTime), .9, .5 * cos(iTime));
}

// 模型
struct Wall {
    vec3 position;
    vec3 normal;
    vec4 albedo;    
} left, right, bottom, top, back, front;

void updateModel() {
    left.position = vec3(1., 0., 0.);
    left.normal = 
}

// 追踪
struct Ray {
    vec3 origin;                                // 射线起点
    vec3 direction;                             // 射线方向
};

struct AABB {
    vec3 min_point;                             // 最小端点
    vec3 max_point;                             // 最大端点
};

vec3 rayAt(const Ray ray, const float t) {

    // 以 origin 为起点 direction 为方向射出 t 的长度
    return ray.origin + t * ray.direction;
}

bool intersectAABB(const Ray ray, const AABB aabb,
    inout float min_t, inout float max_t) {

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

float intersectBox(const Ray ray, const vec3 size, out vec3 normal) {

    // 求出射线与 Cornell Box 的交点
    float min_t = 0., max_t = INFINITY;
    if (intersectAABB(ray, AABB(-size, size), min_t, max_t)) {
        vec3 point = rayAt(ray, min_t) / size;
        
        // 根据分量大小判断交点在哪个面上
        if (abs(point.x) > abs(point.y) && abs(point.x) > abs(point.z)) {
            normal = vec3(point.x > 0. ? 1. : -1., 0., 0.);
        } else if (abs(point.y) > abs(point.x) && abs(point.y) > abs(point.z)) {
            normal = vec3(0., point.y > 0. ? 1. : -1., 0.);
        } else normal = vec3(0., 0., point.z >  0. ? 1. : -1.);
        
        // 返回距离射线出发点更近的那个交点
        return min_t;
    }
    return INFINITY;
}

float intersectPlane(const Ray ray, const vec3 point, const vec3 normal) {

    // dot(orig + t * dir, norm) = dot(point - orig, norm) 
    float denom = dot(ray.direction, normal);
    float t = dot(point - ray.origin, normal) / denom;
    return t > 0. ? t : INFINITY;
}

float intersectLight(const Ray ray) {

    // 求出射线与光源所在平面的交点
    float t = intersectPlane(ray, light.position, light.normal);
    vec3 point = rayAt(ray, t);
    
    // 如果在光源范围内，则与之相交
    if (all(lessThan(abs(light.position - point).xz, vec2(light.size * .5)))) {
        return t;
    }
    return INFINITY;
}

void intersectLeft(const Ray ray, inout float min_t,
    out vec3 point, out vec3 normal, out vec4 albedo) {
    vec3 position = vec3(-1., 0., 0.);
    vec3 leftnorm = vec3(1., 0., 0.);
    vec4 leftAlbedo = 
    float t = intersectPlane(ray, position, leftnorm);
    if (t < min_t) {
        vec3 leftpoint = rayAt(ray, t);
        if (all(lessThanEqual(leftpoint.yz, vec2(1.)))
            && all(greaterThanEqual(leftpoint.yz, vec2(-1.)))) {
            min_t = t;
            point = leftpoint;
            normal = leftnorm;
            
            
            albedo = 
        }
    }
}

float intersect(const Ray ray, inout vec3 point,
    inout vec3 normal, out vec4 albedo) {
    float min_t = INFINITY;
    albedo = vec4(0.);

    float t = intersectLight(ray);
    if (t < min_t) {
        min_t = t;
        point = rayAt(ray, t);
    }

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

    // 为当前片元更新灯光、随机数据和模型
    updateLight();
    updateRandom();
    updateModel();

    // 将 uv 移动至屏幕中间并处理长宽比
    vec2 uv = fragCoord / iResolution.xy - .5;
    float aspect = iResolution.x / iResolution.y;
    if (aspect < 1.) uv.y /= aspect; else uv.x *= aspect;

    // 由此片元发出若干条射线进行路径追踪
    vec3 color = 0.;
    for (int i=0; i<NUM_SAMPLES; i++) {
        vec2 random_result = GetRandom();

        // 将每次追踪采样得到的颜色值累加
        color += sampleColor({

            // 从摄像机出发
            CAMERA_POSITION,

            // 随机方向投射
            normalize(vec3(
                uv+random_result.x*dFdx(uv) + random_result.y*dFdy(uv),
                -1.
            ));
        });
    }

    // 算出平均颜色并做 Gamma 运算得到片元结果
    fragColor = vec4(pow(color/float(NUM_SAMPLES), vec3(1./2.2)), 1.);
}