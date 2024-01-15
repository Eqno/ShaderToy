// 数学
#define PI 3.141592653589                       // 常量 π
#define INFINITY 9999999.                       // 无限大
#define RANDOM_SEED 0u                          // 初始随机种子
#define EPS 1e-5                                // 误差 EPS

struct Random {
    uint seed;                                  // 随机种子
    uint index;                                 // 随机索引
} random;

void updateRandom() {

    // 随机参数 rand seed 和 index
    random = Random(
        RANDOM_SEED,
        uint(dot(gl_FragCoord.xy, vec2(1., 4096.)))
    );
}

vec2 encryptTEA(const uint index, inout uint seed) {

    // 对输入的 arg 进行 TEA 加密
    uvec4 key = uvec4(0xa341316c, 0xc8013ea4, 0xad90777d, 0x7e95761e);
    uint v0 = index, v1 = seed++, sum = 0u, delta = 0x9e3779b9u;
    for (int i=0; i<32; i++) {
        sum += delta;
        v0 += ((v1<<4)+key[0]) ^ (v1+sum) ^ ((v1>>5)+key[1]);
        v1 += ((v0<<4)+key[2]) ^ (v0+sum) ^ ((v0>>5)+key[3]);
    }
    return vec2(v0, v1) / vec2(0xffffffffu);
}

vec2 getRandom() {

    // fract 取小数部分获得随机数
    return fract(encryptTEA(random.index, random.seed));
}

vec2 sampleDisk(const vec2 uv) {

    // 圆碟上采样，uv.x 表示角度，uv.y 表示半径
    float theta = 2. * PI * uv.x;
    float radius = sqrt(uv.y);
    return vec2(cos(theta), sin(theta)) * radius;
}

vec3 sampleCosHemisphere(const vec2 uv) {

    // 半球上采样，圆碟采样的结果，映射到半球上
    vec2 disk = sampleDisk(uv);
    return vec3(disk.x, sqrt(max(0., 1.-dot(disk, disk))), disk.y);
}

mat3 constructONBFrisvad(const vec3 normal) {
	mat3 ret;
	ret[1] = normal;
	if (normal.z < -1. + EPS) {
		ret[0] = vec3(0., -1., 0.);
		ret[2] = vec3(-1., 0., 0.);
	} else {
		float a = 1. / (1. + normal.z);
		float b = -normal.x * normal.y * a;
		ret[0] = vec3(1.0 - normal.x * normal.x * a, b, -normal.x);
		ret[2] = vec3(b, 1.0 - normal.y * normal.y * a, -normal.y);
	}
	return ret;
}

vec3 rotateX(const vec3 p, const float a) {
    float c = cos(a), s = sin(a);
    return vec3(p.x, c*p.y-s*p.z, s*p.y+c*p.z);
}
vec3 rotateY(const vec3 p, const float a) {
    float c = cos(a), s = sin(a);
    return vec3(c*p.x+s*p.z, p.y, -s*p.x+c*p.z);
}
vec3 rotateZ(const vec3 p, const float a) {
    float c = cos(a), s = sin(a);
    return vec3(c*p.x-s*p.y, s*p.x+c*p.y, p.z);
}
vec3 rotate(vec3 point, const vec3 angles) {
    return rotateZ(rotateY(rotateX(point,
        angles.x), angles.y), angles.z);
}

// 光照
struct Material {
    vec3 color;                                 // 表面颜色
    float roughness;                            // 面粗糙度
    vec3 fresnel;                               // 菲涅尔数
    float metallic;                             // 面金属度
    float intensity;                            // 墙体亮度
};
struct Light {
    vec2 size;                                  // 光源大小
    float area;                                 // 光源面积
    vec3 normal;                                // 光源朝向
    vec3 position;                              // 光源位置

    Material mat;
} light;

void updateLight() {
    light.size = vec2(.25, .25);
    light.normal = vec3(0., -1., 0.);
    light.position = vec3(.5 * sin(iTime), .9, .5 * cos(iTime));

    light.area = light.size.x * light.size.y;
    light.mat.color = vec3(1., 1., 1.);
    light.mat.intensity = 2.;
}

// 模型
struct Wall {
    vec2 size;                                  // 墙体大小
    vec3 normal;                                // 墙体法线
    vec3 position;                              // 墙体位置

    Material mat;                               // 墙体材质
} left, right, bottom, top, back;

void updateWalls() {
    left.size = vec2(1., 1.);
    left.normal = vec3(1., 0., 0.);
    left.position = vec3(-1., 0., 0.);

    left.mat.color = vec3(.9, .1, .1);
    left.mat.roughness = .6;
    left.mat.fresnel = vec3(.03, .03, .03);
    left.mat.metallic = .1;
    left.mat.intensity = 0.;

    right.size = vec2(1., 1.);
    right.normal = vec3(-1., 0., 0.);
    right.position = vec3(1, 0., 0.);

    right.mat.color = vec3(.1, .9, .1);
    right.mat.roughness = 1.;
    right.mat.fresnel = vec3(.93, .93, .93);
    right.mat.metallic = 1.;
    right.mat.intensity = 0.;

    bottom.size = vec2(1., 1.);
    bottom.normal = vec3(0., 1., 0.);
    bottom.position = vec3(0., -1., 0.);

    bottom.mat.color = vec3(.7, .7, .7);
    bottom.mat.roughness = 1.;
    bottom.mat.fresnel = vec3(.05, .05, .05);
    bottom.mat.metallic = .1;
    bottom.mat.intensity = 0.;

    top.size = vec2(1., 1.);
    top.normal = vec3(0., -1., 0.);
    top.position = vec3(0., 1., 0.);

    top.mat.color = vec3(.7, .7, .7);
    top.mat.roughness = 1.;
    top.mat.fresnel = vec3(.05, .05, .05);
    top.mat.metallic = .1;
    top.mat.intensity = 0.;

    back.size = vec2(1., 1.);
    back.normal = vec3(0., 0., 1.);
    back.position = vec3(0., 0., -1.);

    back.mat.color = vec3(.7, .7, .7);
    back.mat.roughness = 1.;
    back.mat.fresnel = vec3(.02, .02, .02);
    back.mat.metallic = .1;
    back.mat.intensity = 0.;
}

struct Box {
    vec3 size;
    vec3 rotation;
    vec3 position;

    Material mat;
} big, small;

void updateBoxes() {
    big.size = vec3(.25, .5, .25);
    big.rotation = vec3(0., 1.2, 0.);
    big.position = vec3(-.25, -.5, -.35);

    big.mat.color = vec3(.7, .7, .7);
    big.mat.roughness = 1.;
    big.mat.fresnel = vec3(.03, .03, .03);
    big.mat.metallic = .1;
    big.mat.intensity = 0.;

    small.size = vec3(.25, .25, .25);
    small.rotation = vec3(0., 0., 0.);
    small.position = vec3(.5, -.75, .35);

    small.mat.color = vec3(.7, .7, .7);
    small.mat.roughness = 1.;
    small.mat.fresnel = vec3(.02, .02, .02);
    small.mat.metallic = .1;
    small.mat.intensity = 0.;
}

struct Sphere {
    float radius;
    vec3 position;

    Material mat;
} sphere;

void updateSphere() {
    sphere.radius = .3;
    sphere.position = vec3(-.35, -.75, .5);

    sphere.mat.color = vec3(.7, .7, .7);
    sphere.mat.roughness = .3;
    sphere.mat.fresnel = vec3(1., .71, .29);
    sphere.mat.metallic = .9;
    sphere.mat.intensity = 0.;
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

float _intersectBox(const Ray ray, const vec3 size, out vec3 normal) {

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

#define intersectBox(box) \
{ \
    vec3 surface_normal; \
    Ray ray_transformed = Ray( \
        rotate(ray.origin-box.position, -box.rotation), \
        rotate(ray.direction, -box.rotation) \
    ); \
    float t = _intersectBox(ray_transformed, box.size, surface_normal); \
    if (t < min_t) { \
        min_t = t; \
        point = rayAt(ray, t); \
        normal = rotate(surface_normal, box.rotation); \
        material = box.mat; \
    } \
}

float _intersectPlane(const Ray ray, const vec3 point, const vec3 normal) {

    // dot(orig + t * dir, norm) = dot(point - orig, norm)
    float denom = dot(ray.direction, normal);
    float t = dot(point - ray.origin, normal) / denom;
    return t > 0. ? t : INFINITY;
}

#define intersectPlane(wall, coord) \
{ \
    float t = _intersectPlane(ray, wall.position, wall.normal); \
    if (t < min_t) { \
        vec3 p = rayAt(ray, t); \
        if (all(lessThan(abs(wall.position-p).coord, wall.size))) { \
            min_t = t; \
            point = p; \
            normal = wall.normal; \
            material = wall.mat; \
        } \
    } \
}

float _intersectSphere(const Ray ray, const vec3 center, const float radius) {
    
    vec3 line = center - ray.origin;
    float proj = dot(ray.direction, line);

    float line_square = dot(line, line);
    float radius_square = radius * radius;
    
    if (line_square > radius_square && proj < 0.) {
        return INFINITY;
    }

    float delta_square = line_square - proj * proj;
    if (delta_square > radius_square) {
        return INFINITY;
    }

    float delta = sqrt(radius_square - delta_square);
    if (line_square > radius_square) {
        return proj - delta;
    } else return proj + delta;
}

#define intersectSphere(sphere) \
{ \
    float t = _intersectSphere(ray, sphere.position, sphere.radius); \
    if (t < min_t) { \
        vec3 p = rayAt(ray, t); \
        min_t = t; \
        point = p; \
        normal = normalize(p - sphere.position); \
        material = sphere.mat; \
    } \
}

float intersect(const Ray ray, inout vec3 point,
    inout vec3 normal, out Material material) {

    // 初始化 min_t 为无穷，albedo 为全 0
    float min_t = INFINITY;

    // 求光源与场景中箱子的交点
    intersectBox(big);
    intersectBox(small);
    intersectSphere(sphere);

    // 求光源与墙体和灯光的交点
    intersectPlane(light, xz);
    intersectPlane(left, yz);
    intersectPlane(right, yz);
    intersectPlane(bottom, xz);
    intersectPlane(top, xz);
    intersectPlane(back, xy);

    // 归一化交点处的法向量并返回 t
    normal = normalize(normal);
    return min_t;
}

bool testVisibility(const vec3 p1, const vec3 p2) {

    // 构造 p1 至 p2 的射线
    Ray ray = Ray(p1, normalize(p2 - p1));

    // 使 origin 略微前移防止 p1 干扰
    ray.origin += EPS * ray.direction;
    
    // 求射线与场景的交点
    vec3 normal, point;
    Material material;
    float t = intersect(ray, point, normal, material);

    // 如果 t 大于 p1 与 p2 的距离则 p1、p2 互相可见
    return t > distance(p1, p2) - 2. * EPS;
}

// 材质
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
    vec3 fresnel = mix(mat.fresnel, mat.color, mat.metallic);

    float D = DistributionGGX(dnh, mat.roughness);
    vec3 F = FresnelSchlick(max(dvh, 0.), fresnel);
    float G = GeometrySmith(dnv, dnl, mat.roughness);

    vec3 specular = D * F * G / max(4. * saturate(dnv) * saturate(dnl), EPS);
    vec3 kd = mix(vec3(1.) - F, vec3(0.), mat.metallic);
    return kd * mat.color / PI + specular;
}

// 渲染
#define CAMERA_POSITION vec3(0., 0., 3.125)
#define NUM_SAMPLES 16
#define NUM_BOUNCES 2

vec3 sampleLight(const vec2 offset) {

    // 返回光源上偏离为 offset 的点
    return light.position + vec3(
        (offset.x-.5)*light.size.x,
        0.,
        (offset.y-.5)*light.size.y
    );
}

void directLight(const Ray ray, const vec3 ejectionAtten, const vec3 position,
    const vec3 normal, const Material material, inout vec3 contribution) {

    // 对光源上的点进行随机采样
    vec3 pos_on_light = sampleLight(getRandom());

    // 求出交点和光源点连线的方向
    vec3 direction = pos_on_light - position;
    float product = dot(direction, direction);
    direction /= sqrt(product);

    // 求 BRDF 项
    float brdf_pdf = 1. / PI;
    vec3 brdf = BRDF(normal, direction, -ray.direction, material) / brdf_pdf;

    // 光线出、入因夹角而衰减
    float angularAtten = max(0., dot(normal, direction))
        * max(0., -dot(direction, light.normal)) / product;

    // 如果光亮没有衰减至 0
    if (angularAtten > 0.) {

        // 如果片元处在阴影中，则舍弃直接光照的贡献
        if (testVisibility(position, pos_on_light)) {

            // 渲染方程：弹射衰减 * 角度衰减 * 入射光强 * brdf
            vec3 intensity = light.mat.color * light.mat.intensity;
            contribution += ejectionAtten * angularAtten * intensity * brdf;
        }
    }
}

bool indirectLight(inout Ray ray, inout vec3 ejectionAtten, inout vec3 position,
    inout vec3 normal, inout Material material, inout vec3 contribution) {

    // 利用蒙特卡洛积分确定弹射的射线
    mat3 onb = constructONBFrisvad(normal);
    vec3 direction_next = normalize(onb * sampleCosHemisphere(getRandom()));

    // 求 BRDF 项
    float brdf_pdf = 1. / PI;
    vec3 brdf = BRDF(normal, direction_next, -ray.direction, material) / brdf_pdf;

    // 使 origin 略微前移防止干扰
    Ray ray_next = Ray(position, direction_next); 
    ray_next.origin += ray_next.direction * EPS;

    // 求射线弹射之后和场景的交点
    vec3 position_next, normal_next;
    Material material_next;
    float t = intersect(ray_next, position_next, normal_next, material_next);

    // 如果没有交点，则终止循环
    if (t == INFINITY) return false;

    // 如果交点在光源上，
    if (material_next.intensity > 0.) {

        // 光线出、入因夹角而衰减
        float angularAtten = max(-INFINITY, dot(normal, ray_next.direction))
            * max(-INFINITY, -dot(ray_next.direction, normal_next)) / (t * t);

        // 如果光亮没有衰减至 0
        if (angularAtten > 0.) {

            // 渲染方程：弹射衰减 * 角度衰减 * 入射光强 * brdf
            vec3 intensity = light.mat.color * light.mat.intensity;
            contribution += ejectionAtten * angularAtten * intensity * brdf;
        }
        
        // 停止弹射
        return false;
    }

    // 如果交点不在光源上更新数据
    ray = ray_next;
    position = position_next;
    normal = normal_next;
    material = material_next;

    // 更新弹射造成的光衰减
    ejectionAtten *= brdf;

    // 继续弹射
    return true;
}

vec3 sampleColor(Ray ray) {

    // 初始化光照贡献和弹射衰减
    vec3 contribution = vec3(0.);
    vec3 ejectionAtten = vec3(1.);

    // 求出射线与整个场景的第一个交点
    vec3 position, normal;
    Material material;
    float t = intersect(ray, position, normal, material);

    // 如果不存在交点，返回黑色
    if (t == INFINITY) return vec3(0.);

    // 如果位于光源上，返回光源颜色
    if (material.intensity > 0.) {
        return material.color * material.intensity;
    }

    // 从交点出发，进行射线弹射
    for (int i=0; i<NUM_BOUNCES; i++) {

        // 累加直接光照的贡献值
        directLight(ray, ejectionAtten, position,
            normal, material, contribution);

        // 累加间接光照的贡献值
        if (indirectLight(ray, ejectionAtten, position,
            normal, material, contribution) == false) break;
    }
    return contribution;
}

// 此程序（Pixel Shader）的入口函数
void mainImage(out vec4 fragColor, in vec2 fragCoord) {

    // 为当前片元更新灯光、随机数据和模型
    updateLight();
    updateRandom();
    updateWalls();
    updateBoxes();
    updateSphere();

    // 将 uv 移动至屏幕中间并处理长宽比
    vec2 uv = fragCoord / iResolution.xy - .5;
    float aspect = iResolution.x / iResolution.y;
    if (aspect < 1.) uv.y /= aspect; else uv.x *= aspect;

    // 由此片元发出若干条射线进行路径追踪
    vec3 color = vec3(0.);
    for (int i=0; i<NUM_SAMPLES; i++) {
        vec2 random_result = getRandom();

        // 将每次追踪采样得到的颜色值累加
        color += sampleColor(Ray(

            // 从摄像机出发
            CAMERA_POSITION,

            // 随机方向投射
            normalize(vec3(
                uv+random_result.x*dFdx(uv) + random_result.y*dFdy(uv),
                -1.
            ))
        ));
    }
    
    // 算出平均颜色并做 Gamma 运算得到片元结果
    fragColor = vec4(pow(color/float(NUM_SAMPLES), vec3(1./2.2)), 1.);
}