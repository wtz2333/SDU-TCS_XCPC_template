#include <bits/stdc++.h>
using namespace std;
using ll = long long;
using i128 = __int128;
constexpr double eps = 1e-7;
constexpr double PI = acos(-1);
constexpr double inf = 1e9;
struct Point { double x, y; };        // 点
using Vec = Point;                    // 向量
struct Line { Point P; Vec v; };      // 直线（点向式），射线时为A->B
struct Seg { Point A, B; };           // 线段（存两个端点）
struct Circle { Point O; double r; }; // 圆（存圆心和半径）
using Points = std::vector<Point>;
using ConvexHull = std::vector<Point>;
const Point O = {0, 0};                        // 原点
const Line Ox = {O, {1, 0}}, Oy = {O, {0, 1}}; // 坐标轴

bool eq(double a, double b) { return abs(a - b) < eps; } // ==
bool gt(double a, double b) { return a - b > eps; }      // >
bool lt(double a, double b) { return a - b < -eps; }     // <
bool ge(double a, double b) { return a - b > -eps; }     // >=
bool le(double a, double b) { return a - b < eps; }      // <=
Vec operator + (const Vec &a,const Vec &b){return (Vec){a.x + b.x,a.y + b.y};}
Vec operator - (const Vec &a,const Vec &b){return (Vec){a.x - b.x,a.y - b.y};}
Vec operator * (const Vec &a,const double &b){return (Vec){b * a.x,b * a.y};}
Vec operator * (const double &a,const Vec &b){return (Vec){a * b.x,a * b.y};}
Vec operator / (const Vec &a,const double &b){return (Vec){a.x / b,a.y / b};}
double operator * (const Point &a,const Point &b){return a.x * b.x + a.y * b.y;}// dot // 点乘
double operator ^ (const Point &a,const Point &b){return a.x * b.y - a.y * b.x;}// cross  // 叉乘
bool operator < (const Point& a, const Point& b) {return a.x < b.x || (a.x == b.x && a.y < b.y);}
double len(const Vec &a){return sqrt(a * a);}

ll cross(Point a,Point b){return (ll)a.x * (ll)b.y - (ll)a.y * (ll)b.x;}
ll dot(Point a,Point b){return (ll)a.x * (ll)b.x + (ll)a.y * (ll)b.y;}

double angle(const Vec &a,const Vec &b){return acos(a * b / len(a)/len(b));}

double Polar_angle(Vec &v){return atan2(v.y,v.x);}

int sgn(double x){
    if(abs(x) < eps)
        return 0;
    if(x < 0)
        return -1;
    return 1;
}

Vec r90a(Vec v) { return {-v.y, v.x}; }  // 逆时针旋转90度的向量
Vec r90c(Vec v) { return {v.y, -v.x}; }  // 顺时针旋转90度的向量

// 两向量的夹角余弦
// DEPENDS len, V*V
double cos_t(Vec u, Vec v) { return u * v / len(u) / len(v); }

// 归一化向量（与原向量方向相同的单位向量）
// DEPENDS len
Vec norm(Vec v) { return {v.x / len(v), v.y / len(v)}; }

// 与原向量平行且横坐标大于等于0的单位向量
// DEPENDS d*V, len
Vec pnorm(Vec v) { return (v.x < 0 ? -1 : 1) / len(v) * v; }

// 线段的方向向量
// DEPENDS V-V
// NOTE 直线的方向向量直接访问属性v
Vec dvec(Seg l) { return l.B - l.A; }
//----------------------------------------------------------//
Line line(Point A, Point B) { return {A, B - A}; }

// 斜截式直线
Line line(double k, double b) { return {{0, b}, {1, k}}; }

// 点斜式直线
Line line(Point P, double k) { return {P, {1, k}}; }

// 线段所在直线
// DEPENDS V-V
Line line(Seg l) { return {l.A, l.B - l.A}; }

// 给定直线的横坐标求纵坐标
// NOTE 请确保直线不与y轴平行
double at_x(Line l, double x) { return l.P.y + (x - l.P.x) * l.v.y / l.v.x; }

// 给定直线的纵坐标求横坐标
// NOTE 请确保直线不与x轴平行
double at_y(Line l, double y) { return l.P.x - (y + l.P.y) * l.v.x / l.v.y; }

// 点到直线的垂足
// DEPENDS V-V, V*V, d*V
Point pedal(Point P, Line l) { return l.P - (l.P - P) * l.v / (l.v * l.v) * l.v; }

// 过某点作直线的垂线
// DEPENDS r90c
Line perp(Line l, Point P) { return {P, r90c(l.v)}; }

// 角平分线
// DEPENDS V+V, len, norm
Line bisec(Point P, Vec u, Vec v) { return {P, norm(u) + norm(v)}; }
 
//seg------------------------------------------//

// 线段的方向向量
// DEPENDS V-V
// NOTE 直线的方向向量直接访问属性v
//Vec dvec(Seg l) { return l.B - l.A; }

// 线段中点
Point midp(Seg l) { return {(l.A.x + l.B.x) / 2, (l.A.y + l.B.y) / 2}; }

// 线段中垂线
// DEPENDS r90c, V-V, midp
Line perp(Seg l) { return {midp(l), r90c(l.B - l.A)}; }
//----------------------------------------------------------------//
// 向量是否互相垂直
// DEPENDS eq, V*V
bool verti(Vec u, Vec v) { return eq(u * v, 0); }

// 向量是否互相平行
// DEPENDS eq, V^V
bool paral(Vec u, Vec v) { return eq(u ^ v, 0); }

// 向量是否与x轴平行
// DEPENDS eq
bool paral_x(Vec v) { return eq(v.y, 0); }

// 向量是否与y轴平行
// DEPENDS eq
bool paral_y(Vec v) { return eq(v.x, 0); }

// 点是否在直线上
// DEPENDS eq
bool on(Point P, Line l) { return eq((P.x - l.P.x) * l.v.y, (P.y - l.P.y) * l.v.x); }


// 点是否在射线上
// DEPENDS eq
bool on_ray(Point P, Line l) { return on(P,l) && ((P - l.P) * l.v) >= 0; }

// 点是否在线段上
// DEPENDS eq, len, V-V
bool on(Point P, Seg l) { return eq(len(P - l.A) + len(P - l.B), len(l.A - l.B)); }

// 两个点是否重合
// DEPENDS eq
bool operator==(Point A, Point B) { return eq(A.x, B.x) && eq(A.y, B.y); }

// 两条直线是否重合
// DEPENDS eq, on(L)
bool operator==(Line a, Line b) { return on(a.P, b) && on(a.P + a.v, b); }

// 两条线段是否重合
// DEPENDS eq, P==P
bool operator==(Seg a, Seg b) { return (a.A == b.A && a.B == b.B) || (a.A == b.B && a.B == b.A); }

// 以横坐标为第一关键词、纵坐标为第二关键词比较两个点
// DEPENDS eq, lt
//bool operator<(Point A, Point B) { return lt(A.x, B.x) || (eq(A.x, B.x) && lt(A.y, B.y)); }

// 直线与圆是否相切
// DEPENDS eq, V^V, len
bool tangency(Line l, Circle C) { return eq(abs((C.O ^ l.v) - (l.P ^ l.v)), C.r * len(l.v)); }

// 圆与圆是否相切
// DEPENDS eq, V-V, len
bool tangency(Circle C1, Circle C2) { return eq(len(C1.O - C2.O), C1.r + C2.r); }
//--------------------------------------------------//
// 两点间的距离
// DEPENDS len, V-V
double dis(Point A, Point B) { return len(A - B); }

// 点到直线的距离
// DEPENDS V^V, len
double dis(Point P, Line l) { return abs((P ^ l.v) - (l.P ^ l.v)) / len(l.v); }

// 点到线段的距离
double dis(Point P,Seg l) {
    if(((P - l.A) * (l.B - l.A)) < 0 || ((P - l.B) * (l.A - l.B)) < 0){
        return min(dis(P,l.A),dis(P,l.B));
    }else {
        Line ll = line(l);
        return dis(P,ll);
    }
}
// 平行直线间的距离
// DEPENDS d*V, V^V, len, pnorm
// NOTE 请确保两直线是平行的
double dis(Line a, Line b) { return abs((a.P ^ pnorm(a.v)) - (b.P ^ pnorm(b.v))); }


// 平移
// DEPENDS V+V
Line operator+(Line l, Vec v) { return {l.P + v, l.v}; }
Seg operator+(Seg l, Vec v) { return {l.A + v, l.B + v}; }


// 旋转 逆时针
// DEPENDS V+V, V-V
Point rotate(Point P, double rad) { return {cos(rad) * P.x - sin(rad) * P.y, sin(rad) * P.x + cos(rad) * P.y}; }
Point rotate(Point P, double rad, Point C) { return C + rotate(P - C, rad); }                     // DEPENDS ^1
Line rotate(Line l, double rad, Point C = O) { return {rotate(l.P, rad, C), rotate(l.v, rad)}; }  // DEPENDS ^1, ^2
Seg rotate(Seg l, double rad, Point C = O) { return {rotate(l.A, rad, C), rotate(l.B, rad, C)}; } // DEPENDS ^1, ^2

// 直线与直线交点
// DEPENDS eq, d*V, V*V, V+V, V^V
Points inter(Line a, Line b){
    double c = a.v ^ b.v;
    if (eq(c, 0)) {return {};}
    Vec v = 1 / c * Vec{a.P ^ (a.P + a.v), b.P ^ (b.P + b.v)};
    return {{v * Vec{-b.v.x, a.v.x}, v * Vec{-b.v.y, a.v.y}}};
}

// 线段与线段是否相交
int cross_seg(Seg A,Seg B){
    Point a = A.A,b = A.B,c = B.A,d = B.B;
    ll c1 = (b - a) ^ (c - a),c2 = (b - a) ^ (d - a);
    ll d1 = (d - c) ^ (a - c),d2 = (d - c) ^ (b - c);
    if(sgn(c1) * sgn(c2) < 0 && sgn(d1) * sgn(d2) < 0) {
        return 1;// 严格交
    }
    if(c1 == 0 && sgn((a - c) * (b - c)) <= 0) return 0;
    if(c2 == 0 && sgn((a - d) * (b - d)) <= 0) return 0;
    if(d1 == 0 && sgn((c - a) * (d - a)) <= 0) return 0;
    if(d2 == 0 && sgn((c - b) * (d - b)) <= 0) return 0; // 端点交
    return -1;
}

// 直线与线段相交=>直线与直线相交 + 点是否在线段上
// bool cross_line_seg(Line A,Seg B){
//     Line BB = {B.A,B.B};
//     Points tmp = inter(A,BB);
//     if(tmp.size() == 0)return false;
//     return on(tmp[0],B);
// }
bool cross_line_seg(Line A,Seg B){
    if(abs(A.v ^ (B.A - B.B)) < eps)return false;// 平行
    Vec v1 = B.A - A.P,v2 = B.B - A.P;
    if((v2 ^ v1) < 0){
        swap(v1,v2);
    }else if(abs(v2 ^ v1) < eps){
        if((v1 * v2) <= 0)return true;
        else return false;
    }// 保证 v2 在 v1 下面
    int d1 = sgn(A.v ^ v1);
    int d2 = sgn(A.v ^ v2);
    if(d1 * d2 <= 0)return true;
    return false;
}

// 射线与射线交 
// bool cross_ray_ray(Line A,Line B){  
//     Points tmp = inter(A,B);
//     if(tmp.size() == 0)return false;//注意重合
//     int d1 = sgn((tmp[0] - A.P) * A.v);
//     int d2 = sgn((tmp[0] - B.P) * B.v);
//     return d1 >= 0 && d2 >= 0;
// }

int cross_ray_ray(Line A,Line B) {
    if(fabs(A.v ^ B.v) < eps) {
       if(fabs((A.P - B.P) ^ A.v) < eps) {
           if(sgn(A.v * (B.P - A.P)) < 0 && sgn(B.v * (A.P - B.P)) < 0) return -1;
           else return 0;
       }else return -1;
    }
    Vec v = B.P - A.P;
    double c1 = v ^ A.v;
    double c2 = v ^ B.v;
    double c = A.v ^ B.v;
    if(sgn(c1) * sgn(c) >= 0 && sgn(c2) * sgn(c) >= 0) return 1; // 交
    return -1;
}
// 射线与线段交
// bool cross_ray_seg(Line A,Seg B){
//     Line BB = {B.A,B.B};
//     Points tmp = inter(A,BB);
//     if(tmp.size() == 0)return false;//注意重合
//     int d = sgn((tmp[0] - A.P) * A.v);
//     return on(tmp[0],B) && d >= 0;
// }
bool cross_ray_seg(Line A,Seg B){
    if(abs(A.v ^ (B.A - B.B)) < eps)return false;// 平行
    Vec v1 = B.A - A.P,v2 = B.B - A.P;
    if((v2 ^ v1) < 0){
        swap(v1,v2);
    }else if(abs(v2 ^ v1) < eps){
        if((v1 * v2) <= 0)return true;
        else return false;
    }// 保证 v2 在 v1 下面
    int d1 = sgn(A.v ^ v1);
    int d2 = sgn(A.v ^ v2);
    if(d1 >= 0 && d2 <= 0)return true;
    return false;
}

// 射线与直线交 
// bool cross_ray_line(Line A,Line B) { // A 为射线
//     Points tmp = inter(A,B);
//     if(tmp.size() == 0)return false;
//     int d = sgn((tmp[0] - A.P) * A.v);
//     return d >= 0;
// }
int cross_ray_line(Line A,Line B) { 
    Line B1 = B;
    Line B2 = {B.P,O - B.v};
    int p1 = cross_ray_ray(A,B1);
    int p2 = cross_ray_ray(A,B2);
    if(p1 == 0 || p2 == 0) return 0; // 重合
    else if(p1 == 1 || p2 == 1) return 1; // 交
    return -1;// 不交
}
// 直线与圆交点
// DEPENDS eq, gt, V+V, V-V, V*V, d*V, len, pedal
std::vector<Point> inter(Line l, Circle C){
    Point P = pedal(C.O, l);
    double h = len(P - C.O);
    if (gt(h, C.r)) return {};
    if (eq(h, C.r)) return {P};
    double d = sqrt(C.r * C.r - h * h);
    Vec vec = d / len(l.v) * l.v;
    return {P + vec, P - vec};
}

// 圆与圆的交点
// DEPENDS eq, gt, V+V, V-V, d*V, len, r90c
std::vector<Point> inter(Circle C1, Circle C2){
    Vec v1 = C2.O - C1.O, v2 = r90c(v1);
    double d = len(v1);
    if (gt(d, C1.r + C2.r) || gt(abs(C1.r - C2.r), d)) return {};
    if (eq(d, C1.r + C2.r) || eq(d, abs(C1.r - C2.r))) return {C1.O + C1.r / d * v1};
    double a = ((C1.r * C1.r - C2.r * C2.r) / d + d) / 2;
    double h = sqrt(C1.r * C1.r - a * a);
    Vec av = a / len(v1) * v1, hv = h / len(v2) * v2;
    return {C1.O + av + hv, C1.O + av - hv};
}



// 三角形的重心
Point barycenter(Point A, Point B, Point C){
    return {(A.x + B.x + C.x) / 3, (A.y + B.y + C.y) / 3};
}

// 三角形的外心
// DEPENDS r90c, V*V, d*V, V-V, V+V
// NOTE 给定圆上三点求圆，要先判断是否三点共线
Point circumcenter(Point A, Point B, Point C){
    double a = A * A, b = B * B, c = C * C;
    double d = 2 * (A.x * (B.y - C.y) + B.x * (C.y - A.y) + C.x * (A.y - B.y));
    return 1 / d * r90c(a * (B - C) + b * (C - A) + c * (A - B));
}

// 三角形的内心
// DEPENDS len, d*V, V-V, V+V
Point incenter(Point A, Point B, Point C){
    double a = len(B - C), b = len(A - C), c = len(A - B);
    double d = a + b + c;
    return 1 / d * (a * A + b * B + c * C);
}

// 三角形的垂心
// DEPENDS V*V, d*V, V-V, V^V, r90c
Point orthocenter(Point A, Point B, Point C){
    double n = B * (A - C), m = A * (B - C);
    double d = (B - C) ^ (A - C);
    return 1 / d * r90c(n * (C - B) - m * (C - A));
}
int Point_in_Circumcircle(Point A,Point B,Point C,Point P) {
    if(((B - A) ^ (C - A)) > 0) swap(B,C);
    vector<vector<ll>> f(4);
    f[0].push_back(A.x);f[0].push_back(A.y);f[0].push_back(A * A);f[0].push_back(1);
    f[1].push_back(B.x);f[1].push_back(B.y);f[1].push_back(B * B);f[1].push_back(1);
    f[2].push_back(C.x);f[2].push_back(C.y);f[2].push_back(C * C);f[2].push_back(1);
    f[3].push_back(P.x);f[3].push_back(P.y);f[3].push_back(P * P);f[3].push_back(1);
    vector<int> a(4);
    iota(a.begin(),a.end(),0);
    i128 tot = 0;
    do {    
        i128 sum = 1;
        int op = 1;
        for(int i = 0;i < 4;i ++) {
            sum *= f[i][a[i]];
        }
        for(int i = 0;i < 4;i ++) {
            for(int j = 0;j < i;j ++) {
                if(a[i] < a[j]) op = -op;
            }
        }
        tot += sum * op;
        op = -op;
    }while(next_permutation(a.begin(),a.end()));
    if(tot == 0)return 0; // on
    else if(tot > 0) return 1; // out
    else return -1;// in
}

// Graham扫描法

// DEPENDS eq, lt, cross, V-V, P<P

// double theta(Point p) { return p == O ? -1 / 0. : atan2(p.y, p.x); } // 求极角
// void psort(Points &ps, Point c = O) { // 极角排序
//     sort(ps.begin(), ps.end(), [&](auto a, auto b) {
//         return lt(theta(a - c), theta(b - c));
//     });
// }

//极角排序
int qua(const Point &P) {
  if(P.x == 0 && P.y == 0) return 0;
  if(P.x >= 0 && P.y >= 0) return 1;
  if(P.x < 0&& P.y >= 0) return 2;
  if(P.x < 0 && P.y < 0) return 3;
  if(P.x >= 0 && P.y < 0) return 4;
  exit(-1); 
}
void psort(Points &ps, Point c = O) { // 极角排序
    sort(ps.begin(), ps.end(), [&](auto p1, auto p2) {
        return qua(p1 - c) < qua(p2 - c) || qua(p1 - c) == qua(p2 - c) && ((p1 - c) ^(p2 - c)) > 0;
    });
}

bool check(Point p, Point q, Point r) { // 检查三个点组成的两个向量的旋转方向是否为逆时针
    return lt(0, (q - p) ^ (r - q));
}
ConvexHull Andrew(Points &ps) {
    if(ps.size() == 1){
        return ps;
    }
    sort(ps.begin(), ps.end());
    std::vector<int> I{0}, used(ps.size());
    for (int i = 1; i < ps.size(); i++){
        //std::cout << ps[i].x << " " <<ps[i].y <<"\n";  
        while (I.size() > 1 && !check(ps[I[I.size() - 2]], ps[I.back()], ps[i]))
            used[I.back()] = 0, I.pop_back();
        used[i] = 1, I.push_back(i);
    }// 下凸壳
    int tmp = I.size();
    for (int i = ps.size() - 2; i >= 0; i--){
        if (used[i])
            continue;
        while (I.size() > tmp && !check(ps[I[I.size() - 2]], ps[I.back()], ps[i]))
            used[I.back()] = 0, I.pop_back();
        used[i] = 1, I.push_back(i);
    }// 上凸壳 特别注意上凸壳最后一段如果是垂直x轴的，则不会完全记录进上凸壳内
    Points H;
    for (int i = 0; i < I.size() - 1; i++)
        H.push_back(ps[I[i]]);
    return H;
}//逆时针
ConvexHull chull(Points &ps){
    psort(ps, *min_element(ps.begin(), ps.end())); // 以最左下角的点为极角排序
    Points H{ps[0]};
    for (int i = 1; i < ps.size(); i++){
        while (H.size() > 1 && !check(H[H.size() - 2], H.back(), ps[i]))
            H.pop_back();
        H.push_back(ps[i]);
    }
    return H;
}
ConvexHull operator+(const ConvexHull &A,const ConvexHull B){
    int n = A.size();
    int m = B.size();
    std::vector<Point> v1(n),v2(m);
    for(int i = 0;i < n;i ++){
        v1[i] = A[(i + 1) % n] - A[i];
    }
    for(int i = 0;i < m;i ++){
        v2[i] = B[(i + 1) % m] - B[i];
    }
    ConvexHull C;
    C.push_back(A[0] + B[0]);
    int p1 = 0,p2 = 0;
    while(p1 < n && p2 < m){   
        C.push_back(C.back() + ((v1[p1] ^ v2[p2]) >= 0 ? v1[p1 ++] : v2[p2 ++]));   
    }// 对上凸壳做闵可夫斯基和时将 >= 改为 <= 并且合并凸包时不需要排序
    while(p1 < n)C.push_back(C.back() + v1[p1 ++]);
    while(p2 < m)C.push_back(C.back() + v2[p2 ++]);
    C = Andrew(C);
    return C;
}// 要求凸包起点必须为左下

void test(Points a,Point b){
    int n = a.size();
    int r = 0;
    for(int l = 0;l < n;l ++){
        auto nxt = [&](int x){
            return (x + 1) % n;
        };
        while(nxt(r) != l && check(a[l],a[nxt(r)],b)){
            // b为轴点
            r = nxt(r);
            
            if(l == r)break;
        }
    }//极角排序 转半平面
}

// 半平面交
int sgn(Point a) {
    return a.y > 0 || (a.y == 0 && a.x > 0) ? 1 : -1;
}
bool pointOnLineLeft(Point p, Line l) {
    return (l.v ^ (p - l.P)) > eps;
}
Point lineIntersection(Line l1, Line l2) {
    return l1.P + l1.v * ((l2.v ^ (l1.P - l2.P)) / (l2.v ^ (O - l1.v)));
}
std::vector<Point> hp(std::vector<Line> lines) {
    std::sort(lines.begin(), lines.end(), [&](auto l1, auto l2) {
        auto d1 = l1.v;
        auto d2 = l2.v;
         
        if (sgn(d1) != sgn(d2)) {
            return sgn(d1) == 1;
        }
         
        return (d1 ^ d2) > 0;
    });
    std::deque<Line> ls;
    std::deque<Point> ps;
    for (auto l : lines) {
        if (ls.empty()) {
            ls.push_back(l);
            continue;
        }
        while (!ps.empty() && !pointOnLineLeft(ps.back(), l)) {
            ps.pop_back();ls.pop_back();
        }
         
        while (!ps.empty() && !pointOnLineLeft(ps[0], l)) {
            ps.pop_front();ls.pop_front();
        }
        if (abs(cross(l.v, ls.back().v)) < eps) {

            if ((l.v * ls.back().v) > eps) {
                //continue;
                if (!pointOnLineLeft(ls.back().P, l)) {
                    assert(ls.size() == 1);
                    ls[0] = l;
                }
                continue;
            }
            return {};
        }
        auto now = inter(ls.back(),l);
        ps.push_back(now[0]);
        // ps.push_back(lineIntersection(ls.back(), l));
        ls.push_back(l);
    }
     
    while (!ps.empty() && !pointOnLineLeft(ps.back(), ls[0])) {
        ps.pop_back();ls.pop_back();
    }
    if (ls.size() <= 2) {
        return {};
    }
    auto now = inter(ls[0],ls.back());
    ps.push_back(now[0]);
    // ps.push_back(lineIntersection(ls[0], ls.back()));
    return std::vector(ps.begin(), ps.end());
}//逆时针平面，使用时，可以给最外面套个值域的正方形框，来减少精度错误

double get_longest(vector<Point> A) {  // 求凸包直径
  int j = 2,n = A.size();
  double mx = 0;
  if(n == 1) return 0;
  if (n < 3) {
    mx = dis(A[0], A[1]);
    return mx;
  }
  auto sqr = [&](Point A,Point B,Point C) {
    return abs((B - A) ^ (C - B));
  };
  for(int i = 0;i < A.size();i ++) {
    while (sqr(A[i], A[(i + 1) % n], A[j]) <=
           sqr(A[i], A[(i + 1) % n], A[(j + 1) % n]))
      j = (j + 1) % n;
    mx = max(mx, max(dis(A[(i + 1) % n], A[j]), dis(A[i], A[j])));
  }
  return mx;
}

bool check_v_in_AB(Vec v,Vec A,Vec B) {// 判断向量v在向量A和B之间
    if(sgn(A ^ B) < 0)swap(A,B);
    return sgn(A ^ v) >= 0 && sgn(v ^ B) >= 0; 
}
bool Point_in_Triangle(Point x,Point A,Point B,Point C) {
    if(check_v_in_AB(x - A,B - A,C - A) && check_v_in_AB(x - B,A - B,C - B) && check_v_in_AB(x - C,A - C,B - C)) {
        return true;
    }return false;
}
bool Point_in_ConvexHull(Point A,ConvexHull &a) {
    int n = a.size();
    if(n < 3) exit(-1);
    int l = 2,r = n - 1;
    int x = -1;
    while(l <= r) {
        int mid = l + r >> 1;
        if(check_v_in_AB(A - a[0],a[1] - a[0],a[mid] - a[0])) {
            r = mid - 1;
            x = mid;
        }else l = mid + 1;
    }
    if(x == -1) return false;
    return Point_in_Triangle(A,a[0],a[x - 1],a[x]);
}