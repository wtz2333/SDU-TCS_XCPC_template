using LL = long long;
using point_t = double; //全局数据类型，可修改为 long long 等

const point_t eps = 1e-8;
const double PI = acosl(-1);
//const double PI = numbers::pi_v<double>;

// 点与向量
template <typename T>
struct point{
    T x, y;

    bool operator==(const point &a) const { return (abs(x - a.x) <= eps && abs(y - a.y) <= eps); }
    bool operator<(const point &a) const{
        if (abs(x - a.x) <= eps)
            return y < a.y - eps;
        return x < a.x - eps;
    }
    bool operator>(const point &a) const { return !(*this < a || *this == a); }
    point operator+(const point &a) const { return {x + a.x, y + a.y}; }
    point operator-(const point &a) const { return {x - a.x, y - a.y}; }
    point operator-() const { return {-x, -y}; }
    point operator*(const T k) const { return {k * x, k * y}; }
    point operator/(const T k) const { return {x / k, y / k}; }
    T operator*(const point &a) const { return x * a.x + y * a.y; } // 点积
    T operator^(const point &a) const { return x * a.y - y * a.x; } // 叉积，注意优先级
    int toleft(const point &a) const
    {
        const auto t = (*this) ^ a;
        return (t > eps) - (t < -eps);
    }                                                             // to-left 测试
    T len2() const { return (*this) * (*this); }                  // 向量长度的平方
    T dis2(const point &a) const { return (a - (*this)).len2(); } // 两点距离的平方

    // 涉及浮点数
    double len() const { return sqrtl(len2()); }                                                                      // 向量长度
    double dis(const point &a) const { return sqrtl(dis2(a)); }                                                       // 两点距离
    double ang(const point &a) const { return acosl(max(-1.0l, min(1.0l, ((*this) * a) / (len() * a.len())))); }      // 向量夹角
    point rot(const double rad) const { return {x * cos(rad) - y * sin(rad), x * sin(rad) + y * cos(rad)}; }          // 逆时针旋转（给定角度）
    point rot(const double cosr, const double sinr) const { return {x * cosr - y * sinr, x * sinr + y * cosr}; } // 逆时针旋转（给定角度的正弦与余弦）
};

using Point = point<point_t>;

// 极角排序
struct argcmp
{
    bool operator()(const Point &a, const Point &b) const
    {
        const auto quad = [](const Point &a)
        {
            if (a.y < -eps)
                return 1;
            if (a.y > eps)
                return 4;
            if (a.x < -eps)
                return 5;
            if (a.x > eps)
                return 3;
            return 2;
        };
        const int qa = quad(a), qb = quad(b);
        if (qa != qb)
            return qa < qb;
        const auto t = a ^ b;
        // if (abs(t)<=eps) return a*a<b*b-eps;  // 不同长度的向量需要分开
        return t > eps;
    }
};



// 直线
template <typename T>
struct line
{
    point<T> p, v; // p 为直线上一点，v 为方向向量

    bool operator==(const line &a) const { return v.toleft(a.v) == 0 && v.toleft(p - a.p) == 0; }
    int toleft(const point<T> &a) const { return v.toleft(a - p); } // to-left 测试
    bool operator<(const line &a) const                             // 半平面交算法定义的排序
    {
        if (abs(v ^ a.v) <= eps && v * a.v >= -eps)
            return toleft(a.p) == -1;
        return argcmp()(v, a.v);
    }

    // 涉及浮点数
    point<T> inter(const line &a) const { return p + v * ((a.v ^ (p - a.p)) / (v ^ a.v)); } // 直线交点
    double dis(const point<T> &a) const { return abs(v ^ (a - p)) / v.len(); }         // 点到直线距离
    point<T> proj(const point<T> &a) const { return p + v * ((v * (a - p)) / (v * v)); }    // 点在直线上的投影
    point<T> symmetry(const point<T> &a) const { return proj(a) * 2 - a;}                    // 点关于直线的对称点 
};

using Line = line<point_t>;

//线段
template <typename T>
struct segment
{
    point<T> a, b;

    bool operator<(const segment &s) const { return make_pair(a, b) < make_pair(s.a, s.b); }

    // 判定性函数建议在整数域使用

    // 判断点是否在线段上
    // -1 点在线段端点 | 0 点不在线段上 | 1 点严格在线段上
    int is_on(const point<T> &p) const
    {
        if (p == a || p == b)
            return -1;
        return (p - a).toleft(p - b) == 0 && (p - a) * (p - b) < -eps;
    }

    // 判断线段直线是否相交
    // -1 直线经过线段端点 | 0 线段和直线不相交 | 1 线段和直线严格相交
    int is_inter(const line<T> &l) const
    {
        if (l.toleft(a) == 0 || l.toleft(b) == 0)
            return -1;
        return l.toleft(a) != l.toleft(b);
    }

    // 判断两线段是否相交
    // -1 在某一线段端点处相交 | 0 两线段不相交 | 1 两线段严格相交
    int is_inter(const segment<T> &s) const
    {
        if (is_on(s.a) || is_on(s.b) || s.is_on(a) || s.is_on(b))
            return -1;
        const line<T> l{a, b - a}, ls{s.a, s.b - s.a};
        return l.toleft(s.a) * l.toleft(s.b) == -1 && ls.toleft(a) * ls.toleft(b) == -1;
    }

    // 点到线段距离
    double dis(const point<T> &p) const
    {
        if ((p - a) * (b - a) < -eps || (p - b) * (a - b) < -eps)
            return min(p.dis(a), p.dis(b));
        const line<T> l{a, b - a};
        return l.dis(p);
    }

    // 两线段间距离
    double dis(const segment<T> &s) const
    {
        if (is_inter(s))
            return 0;
        return min({dis(s.a), dis(s.b), s.dis(a), s.dis(b)});
    }

    // 只求整点交点可以不使用浮点数,避免精度问题,使用前需要先判断是否有交点 
    pair<bool, point<T> > int_inter(const segment &s){

        // 线段转为直线的一般式
        auto seg2line = [&](const segment &s){
            T A = s.a.y - s.b.y;
            T B = s.b.x - s.a.x;
            T C = -A * s.a.x - B * s.a.y;
            return array{A, B, C};
        };

        auto [A1, B1, C1] = seg2line(*this);
        auto [A2, B2, C2] = seg2line(s);
        T dx = C1 * B2 - C2 * B1;
        T dy = A1 * C2 - A2 * C1;
        T d = B1 * A2 - B2 * A1;
        if (d == 0) return {false, {}};
        if (dy % d || dx % d) return {false, {}};
        return {true, {dx / d, dy / d}};
    }

};
using Segment = segment<point_t>;

template <typename T>
struct polygon
{
    vector<point<T>> p; // 以逆时针顺序存储

    size_t nxt(const size_t i) const { return i == p.size() - 1 ? 0 : i + 1; }
    size_t pre(const size_t i) const { return i == 0 ? p.size() - 1 : i - 1; }

    // 回转数
    // 返回值第一项表示点是否在多边形边上
    // 对于狭义多边形，回转数为 0 表示点在多边形外，否则点在多边形内
    pair<bool, int> winding(const point<T> &a) const
    {
        int cnt = 0;
        for (size_t i = 0; i < p.size(); i++)
        {
            const point<T> u = p[i], v = p[nxt(i)];
            if (abs((a - u) ^ (a - v)) <= eps && (a - u) * (a - v) <= eps)
                return {true, 0};
            if (abs(u.y - v.y) <= eps)
                continue;
            const Line uv = {u, v - u};
            if (u.y < v.y - eps && uv.toleft(a) <= 0)
                continue;
            if (u.y > v.y + eps && uv.toleft(a) >= 0)
                continue;
            if (u.y < a.y - eps && v.y >= a.y - eps)
                cnt++;
            if (u.y >= a.y - eps && v.y < a.y - eps)
                cnt--;
        }
        return {false, cnt};
    }

    // 射线法 2表示在多边形内,1表示在多边形上,0表示在多边形外
    int is_in(const point<T> &a){
        int x = 0;
        for(size_t i = 0; i < p.size(); i++){
            segment<T> s = {p[i], p[nxt(i)]};
            if (s.is_on(a)) return 1;
            point<T> p1 = p[i] - a, p2 = p[nxt(i)] - a;
            if(p1.y > p2.y) swap(p1, p2);
            if(p1.y < eps && p2.y > eps && (p1 ^ p2) > eps) x = !x;
        }
        return x ? 2 : 0;
    }

    // 多边形面积的两倍
    // 可用于判断点的存储顺序是顺时针或逆时针
    T area() const
    {
        T sum = 0;
        for (size_t i = 0; i < p.size(); i++)
            sum += p[i] ^ p[nxt(i)];
        return sum;
    }

    // 多边形的周长
    double circ() const
    {
        double sum = 0;
        for (size_t i = 0; i < p.size(); i++)
            sum += p[i].dis(p[nxt(i)]);
        return sum;
    }
};

using Polygon = polygon<point_t>;

using Point = point<point_t>;
struct Circle
{
    Point c;
    double r;

    bool operator==(const Circle &a) const { return c == a.c && abs(r - a.r) <= eps; }
    double circ() const { return 2 * PI * r; } // 周长
    double area() const { return PI * r * r; } // 面积

    // 点与圆的关系
    // -1 圆上 | 0 圆外 | 1 圆内
    int is_in(const Point &p) const
    {
        const double d = p.dis(c);
        return abs(d - r) <= eps ? -1 : d < r - eps;
    }

    // 直线与圆关系
    // 0 相离 | 1 相切 | 2 相交
    int relation(const Line &l) const
    {
        const double d = l.dis(c);
        if (d > r + eps)
            return 0;
        if (abs(d - r) <= eps)
            return 1;
        return 2;
    }

    // 圆与圆关系
    // -1 相同 | 0 相离 | 1 外切 | 2 相交 | 3 内切 | 4 内含
    int relation(const Circle &a) const
    {
        if (*this == a)
            return -1;
        const double d = c.dis(a.c);
        if (d > r + a.r + eps)
            return 0;
        if (abs(d - r - a.r) <= eps)
            return 1;
        if (abs(d - abs(r - a.r)) <= eps)
            return 3;
        if (d < abs(r - a.r) - eps)
            return 4;
        return 2;
    }

    // 直线与圆的交点
    vector<Point> inter(const Line &l) const
    {
        const double d = l.dis(c);
        const Point p = l.proj(c);
        const int t = relation(l);
        if (t == 0)
            return vector<Point>();
        if (t == 1)
            return vector<Point>{p};
        const double k = sqrt(r * r - d * d);
        return vector<Point>{p - (l.v / l.v.len()) * k, p + (l.v / l.v.len()) * k};
    }

    // 圆与圆交点
    vector<Point> inter(const Circle &a) const
    {
        const double d = c.dis(a.c);
        const int t = relation(a);
        if (t == -1 || t == 0 || t == 4)
            return vector<Point>();
        Point e = a.c - c;
        e = e / e.len() * r;
        if (t == 1 || t == 3)
        {
            if (r * r + d * d - a.r * a.r >= -eps)
                return vector<Point>{c + e};
            return vector<Point>{c - e};
        }
        const double costh = (r * r + d * d - a.r * a.r) / (2 * r * d), sinth = sqrt(1 - costh * costh);
        return vector<Point>{c + e.rot(costh, -sinth), c + e.rot(costh, sinth)};
    }

    // 圆与圆交面积
    double inter_area(const Circle &a) const
    {
        const double d = c.dis(a.c);
        const int t = relation(a);
        if (t == -1)
            return area();
        if (t < 2)
            return 0;
        if (t > 2)
            return min(area(), a.area());
        const double costh1 = (r * r + d * d - a.r * a.r) / (2 * r * d), costh2 = (a.r * a.r + d * d - r * r) / (2 * a.r * d);
        const double sinth1 = sqrt(1 - costh1 * costh1), sinth2 = sqrt(1 - costh2 * costh2);
        const double th1 = acos(costh1), th2 = acos(costh2);
        return r * r * (th1 - costh1 * sinth1) + a.r * a.r * (th2 - costh2 * sinth2);
    }

    // 过圆外一点圆的切线
    vector<Line> tangent(const Point &a) const
    {
        const int t = is_in(a);
        if (t == 1)
            return vector<Line>();
        if (t == -1)
        {
            const Point v = {-(a - c).y, (a - c).x};
            return vector<Line>{{a, v}};
        }
        Point e = a - c;
        e = e / e.len() * r;
        const double costh = r / c.dis(a), sinth = sqrt(1 - costh * costh);
        const Point t1 = c + e.rot(costh, -sinth), t2 = c + e.rot(costh, sinth);
        return vector<Line>{{a, t1 - a}, {a, t2 - a}};
    }

    // 两圆的公切线
    vector<Line> tangent(const Circle &a) const
    {
        const int t = relation(a);
        vector<Line> lines;
        if (t == -1 || t == 4)
            return lines;
        if (t == 1 || t == 3)
        {
            const Point p = inter(a)[0], v = {-(a.c - c).y, (a.c - c).x};
            lines.push_back({p, v});
        }
        const double d = c.dis(a.c);
        const Point e = (a.c - c) / (a.c - c).len();
        if (t <= 2)
        {
            const double costh = (r - a.r) / d, sinth = sqrt(1 - costh * costh);
            const Point d1 = e.rot(costh, -sinth), d2 = e.rot(costh, sinth);
            const Point u1 = c + d1 * r, u2 = c + d2 * r, v1 = a.c + d1 * a.r, v2 = a.c + d2 * a.r;
            lines.push_back({u1, v1 - u1});
            lines.push_back({u2, v2 - u2});
        }
        if (t == 0)
        {
            const double costh = (r + a.r) / d, sinth = sqrt(1 - costh * costh);
            const Point d1 = e.rot(costh, -sinth), d2 = e.rot(costh, sinth);
            const Point u1 = c + d1 * r, u2 = c + d2 * r, v1 = a.c - d1 * a.r, v2 = a.c - d2 * a.r;
            lines.push_back({u1, v1 - u1});
            lines.push_back({u2, v2 - u2});
        }
        return lines;
    }

    // 圆的反演
    tuple<int, Circle, Line> inverse(const Line &l) const
    {
        const Circle null_c = {{0.0, 0.0}, 0.0};
        const Line null_l = {{0.0, 0.0}, {0.0, 0.0}};
        if (l.toleft(c) == 0)
            return {2, null_c, l};
        const Point v = l.toleft(c) == 1 ? Point{l.v.y, -l.v.x} : Point{-l.v.y, l.v.x};
        const double d = r * r / l.dis(c);
        const Point p = c + v / v.len() * d;
        return {1, {(c + p) / 2, d / 2}, null_l};
    }

    tuple<int, Circle, Line> inverse(const Circle &a) const
    {
        const Circle null_c = {{0.0, 0.0}, 0.0};
        const Line null_l = {{0.0, 0.0}, {0.0, 0.0}};
        const Point v = a.c - c;
        if (a.is_in(c) == -1)
        {
            const double d = r * r / (a.r + a.r);
            const Point p = c + v / v.len() * d;
            return {2, null_c, {p, {-v.y, v.x}}};
        }
        if (c == a.c)
            return {1, {c, r * r / a.r}, null_l};
        const double d1 = r * r / (c.dis(a.c) - a.r), d2 = r * r / (c.dis(a.c) + a.r);
        const Point p = c + v / v.len() * d1, q = c + v / v.len() * d2;
        return {1, {(p + q) / 2, p.dis(q) / 2}, null_l};
    }
};

Circle get_circle(Point a, Point b, Point c){
    Line l1 = {(a + b) / 2, (a - b).rot(PI / 2)};
    Line l2 = {(a + c) / 2, (a - c).rot(PI / 2)};
    Point p = l1.inter(l2);
    return {p, p.dis(a)};
}
