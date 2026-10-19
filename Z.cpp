using namespace std;
using i64 = long long;
constexpr int P = 998244353;
int norm(int x) {
    if (x < 0) {
        x += P;
    }
    if (x >= P) {
        x -= P;
    }
    return x;
}
template<class T>
T qp(T a, int b) {
    T res = 1;
    for (; b; b /= 2, a *= a) {
        if (b % 2) {
            res *= a;
        }
    }
    return res;
}
struct Z{
    int x;
    Z(): x{} {}
    Z(int x) : x{norm(x)} {}
    Z(i64 x) : x{norm(int(x % P))} {}
    friend std::istream &operator>>(std::istream &is,Z &a) {
        i64 v;
        is >> v;
        a = Z(v);
        return is;
    }
    friend std::ostream &operator<<(std::ostream &os, const Z &a) {
        return os << a.x;
    }
    Z inv() const {
        return qp(Z(x),P - 2);
    }
};
bool operator==(const Z a,const Z b) { return a.x == b.x; }
bool operator!=(const Z a,const Z b) { return a.x != b.x; }
Z operator+(const Z a, const Z b) { return norm(a.x + b.x); }
Z operator-(const Z a, const Z b) { return norm(a.x + P - b.x); }
Z operator-(const Z x) { return x.x ? P - x.x : 0; }
Z operator*(const Z a, const Z b) { return i64(a.x) * b.x % P; }
Z operator/(const Z a, const Z b) { return a * b.inv(); }
Z &operator+=(Z &a, const Z b) { return a = a + b; }
Z &operator-=(Z &a, const Z b) { return a = a - b; }
Z &operator*=(Z &a, const Z b) { return a = a * b; }
Z &operator/=(Z &a, const Z b) { return a = a / b; }