// 数学
#define PI 3.141592653589                       // 常量 π
#define INFINITY 9999999.                       // 无限大
#define RANDOM_SEED 0                           // 初始随机种子
#define EPS 1e-5                                // 误差 EPS

struct Random {
    int seed;                                   // 随机种子
    int flat_index;                             // 随机索引
} random;

void updateRandom() {

    // 随机参数 rand seed 和 index
    random = Random(
        RANDOM_SEED,
        int(dot(gl_FragCoord.xy, vec2(1., 4096.)))
    );
}

void EncryptTEA(inout uvec2 arg) {

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
    EncryptTEA(arg);
    return fract(vec2(arg) / vec2(0xffffffffu));
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
	if(normal.z < -0.999805696) {
		ret[0] = vec3(0.0, -1.0, 0.0);
		ret[2] = vec3(-1.0, 0.0, 0.0);
	}
	else {
		float a = 1.0 / (1.0 + normal.z);
		float b = -normal.x * normal.y * a;
		ret[0] = vec3(1.0 - normal.x * normal.x * a, b, -normal.x);
		ret[2] = vec3(b, 1.0 - normal.y * normal.y * a, -normal.y);
	}
	return ret;
}

mat4 matRotateX(const float a) {
	mat4 ret = mat4(1.);
	ret[1][1] = ret[2][2] = cos(a);
	ret[2][1] = sin(a);
	ret[1][2] = -ret[2][1];
	return ret;
}
mat4 matRotateY(const float a) {
	mat4 ret = mat4(1.);
	ret[0][0] = ret[2][2] = cos(a);
	ret[0][2] = sin(a);
	ret[2][0] = -ret[0][2];
	return ret;
}
mat4 matRotateZ(const float a) {
	mat4 ret = mat4(1.);
	ret[0][0] = ret[1][1] = cos(a);
	ret[1][0] = sin(a);
	ret[0][1] = -ret[1][0];
	return ret;
}
vec3 rotate(const vec4 vec, const vec3 angles) {
    return vec3(matRotateX(angles.x) * matRotateY(angles.y)
        * matRotateZ(angles.z) * vec);
}

// 光照
struct Light {
    vec2 size;                                  // 光源大小
    vec3 normal;                                // 光源朝向
    vec3 position;                              // 光源位置

    vec3 color;                                 // 光源颜色
    float intensity;                            // 光源亮度
    vec4 albedo;                                // 光反射率
    float area;                                 // 光源区域
} light;

void updateLight() {
    light.size = vec2(.25, .25);
    light.normal = vec3(0., -1., 0.);
    light.position = vec3(.5 * sin(iTime), .9, .5 * cos(iTime));

    light.color = vec3(1., 1., 1.);
    light.intensity = 2.;

    light.area = light.size.x * light.size.y;
    light.albedo = vec4(light.color, light.intensity / light.area);
}

// 模型
struct Wall {
    vec2 size;                                  // 墙体大小
    vec3 normal;                                // 墙体法线
    vec3 position;                              // 墙体位置

    vec3 color;                                 // 墙体颜色
    float intensity;                            // 墙体亮度
    vec4 albedo;                                // 光反射率
} left, right, bottom, top, back;

void updateWalls() {
    left.size = vec2(1., 1.);
    left.normal = vec3(1., 0., 0.);
    left.position = vec3(-1., 0., 0.);

    left.color = vec3(.9, .1, .1);
    left.intensity = 0.;
    left.albedo = vec4(left.color, left.intensity);

    right.size = vec2(1., 1.);
    right.normal = vec3(-1., 0., 0.);
    right.position = vec3(1., 0., 0.);

    right.color = vec3(.1, .9, .1);
    right.intensity = 0.;
    right.albedo = vec4(right.color, right.intensity);

    bottom.size = vec2(1., 1.);
    bottom.normal = vec3(0., 1., 0.);
    bottom.position = vec3(0., -1., 0.);

    bottom.color = vec3(.7, .7, .7);
    bottom.intensity = 0.;
    bottom.albedo = vec4(bottom.color, bottom.intensity);

    top.size = vec2(1., 1.);
    top.normal = vec3(0., -1., 0.);
    top.position = vec3(0., 1., 0.);

    top.color = vec3(.7, .7, .7);
    top.intensity = 0.;
    top.albedo = vec4(top.color, top.intensity);

    back.size = vec2(1., 1.);
    back.normal = vec3(0., 0., 1.);
    back.position = vec3(0., 0., -1.);

    back.color = vec3(.7, .7, .7);
    back.intensity = 0.;
    back.albedo = vec4(back.color, back.intensity);

}

struct Box {
    vec3 size;
    vec3 rotation;
    vec3 position;

    vec3 color;
    float intensity;
    vec4 albedo;
} big, small;

void updateBoxes() {
    big.size = vec3(.25, .5, .25);
    big.rotation = vec3(0., .6, 0.);
    big.position = vec3(-.35, -.5, -.35);

    big.color = vec3(.7, .7, .7);
    big.intensity = 0.;
    big.albedo = vec4(big.color, big.intensity);

    small.size = vec3(0.25, 0.25, 0.25);
    small.rotation = vec3(0., 0., 0.);
    small.position = vec3(0.5, -0.75, 0.35);

    small.color = vec3(.7, .7, .7);
    small.intensity = 0.;
    small.albedo = vec4(small.color, small.intensity);
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
        rotate(vec4(ray.origin-box.position, 1.), box.rotation), \
        rotate(vec4(ray.direction, 0.), box.rotation) \
    ); \
    float t = _intersectBox(ray_transformed, box.size, surface_normal); \
    if (t < min_t) { \
        min_t = t; \
        point = rayAt(ray, t); \
        albedo = box.albedo; \
        normal = rotate(vec4(surface_normal, 0.), -box.rotation); \
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
        if (all(lessThan(abs(wall.position-p).coord, vec2(wall.size)))) { \
            min_t = t; \
            point = p; \
            normal = wall.normal; \
            albedo = wall.albedo; \
        } \
    } \
}

float intersect(const Ray ray, inout vec3 point,
    inout vec3 normal, out vec4 albedo) {

    // 初始化 min_t 为无穷，albedo 为全 0
    float min_t = INFINITY;
    albedo = vec4(0.);

    // 求光源与场景中箱子的交点
    intersectBox(big);
    intersectBox(small);

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
    vec4 albedo;
    float t = intersect(ray, point, normal, albedo);

    // 如果 t 大于 p1 与 p2 的距离则 p1、p2 互相可见
    return t > distance(p1, p2) - 2. * EPS;
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

void directLight(const vec3 ejectionAtten, const vec3 position,
    const vec3 normal, const vec4 albedo, inout vec3 contribution) {

    // 求 brdf_pdf 和 brdf
    float brdf_pdf = 1. / PI;
    vec3 brdf = albedo.rgb / PI;

    // 对光源上的点进行随机采样
    vec3 pos_on_light = sampleLight(GetRandom());

    // 求出交点和光源点连线的方向
    vec3 direction = pos_on_light - position;
    float product = dot(direction, direction);
    direction /= sqrt(product);

    // 光线出、入因夹角而衰减
    float angularAtten = max(0., dot(normal, direction))
        * max(0., -dot(direction, light.normal)) / product;

    // 如果光亮没有衰减至 0
    if (angularAtten > 0.) {

        // 求算直接光照 PDF 在所有光照 PDF 中的比重
        float light_pdf = 1. / (light.area * angularAtten);
        float weight = light_pdf / (light_pdf + brdf_pdf);

        // 如果片元处在阴影中，则舍弃直接光照的贡献
        if (testVisibility(position, pos_on_light)) {

            // 渲染方程：弹射衰减 * 入射光强 * 直接比重 (* 角度衰减) / 光照分布。
            vec3 intensity = light.albedo.rgb * light.albedo.a;
            contribution += ejectionAtten * intensity * brdf / light_pdf;
        }
    }
}

bool indirectLight(inout vec3 ejectionAtten, inout vec3 position,
    inout vec3 normal, inout vec4 albedo, inout vec3 contribution) {

    // 求 brdf_pdf 和 brdf
    float brdf_pdf = 1. / PI;
    vec3 brdf = albedo.rgb / PI;

    // 利用蒙特卡洛积分确定弹射的射线
    mat3 onb = constructONBFrisvad(normal);
    vec3 direction_next = normalize(onb * sampleCosHemisphere(GetRandom()));

    // 使 origin 略微前移防止干扰
    Ray ray_next = Ray(position, direction_next); 
    ray_next.origin += ray_next.direction * EPS;

    // 求射线弹射之后和场景的交点
    vec3 position_next, normal_next;
    vec4 albedo_next;
    float t = intersect(ray_next, position_next, normal_next, albedo_next);

    // 如果没有交点，则并终止循环
    if (t == INFINITY) return false;

    // 如果交点在光源上，
    if (albedo_next.a > 0.) {

        // 光线出、入因夹角而衰减
        float angularAtten = max(0., dot(normal, ray_next.direction))
            * max(0., -dot(ray_next.direction, normal_next)) / (t * t);

        // 如果光亮没有衰减至 0
        if (angularAtten > 0.) {

            // 求算间接光照 PDF 在所有光照 PDF 中的比重
            float light_pdf = 1. / (light.area * angularAtten);
            float weight = brdf_pdf / (light_pdf + brdf_pdf);

            // 渲染方程：弹射衰减 * 入射光强 * 直接比重 (* 角度衰减) / 光照分布。
            vec3 intensity = light.albedo.rgb * light.albedo.a;
            contribution += ejectionAtten * intensity  * brdf / light_pdf;
        }
        
        // 停止弹射
        return false;
    }

    // 如果交点不在光源上，更新起点为交点
    ejectionAtten *= brdf / brdf_pdf;
    position = position_next;
    normal = normal_next;
    albedo = albedo_next;

    // 继续弹射
    return true;
}

vec3 sampleColor(const Ray ray) {
    // 初始化光照贡献和弹射衰减
    vec3 contribution = vec3(0.);
    vec3 ejectionAtten = vec3(1.);

    // 求出射线与整个场景的第一个交点
    vec3 position, normal;
    vec4 albedo;
    float t = intersect(ray, position, normal, albedo);

    // 如果不存在交点，返回黑色
    if (t == INFINITY) return vec3(0.);

    // 如果位于光源上，返回光源颜色
    if (albedo.a > 0.) return albedo.rgb * albedo.a;

    // 从交点出发，进行射线弹射
    for (int i=0; i<NUM_BOUNCES; i++) {

        // 累加直接光照的贡献值
        // directLight(ejectionAtten, position,
        //     normal, albedo, contribution);

        // 累加间接光照的贡献值
        if (indirectLight(ejectionAtten, position,
            normal, albedo, contribution) == false) break;
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

    // 将 uv 移动至屏幕中间并处理长宽比
    vec2 uv = fragCoord / iResolution.xy - .5;
    float aspect = iResolution.x / iResolution.y;
    if (aspect < 1.) uv.y /= aspect; else uv.x *= aspect;

    // 由此片元发出若干条射线进行路径追踪
    vec3 color = vec3(0.);
    for (int i=0; i<NUM_SAMPLES; i++) {
        vec2 random_result = GetRandom();

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