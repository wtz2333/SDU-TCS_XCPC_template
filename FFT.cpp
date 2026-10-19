#include<bits/stdc++.h>
using namespace std;
typedef long long ll;
constexpr double PI = std::atan2(0, -1);
std::vector<int> rev;
std::vector<std::complex<double>> roots {0, 1};
void dft(std::vector<std::complex<double>> &a) {
    int n = a.size();
    if (int(rev.size()) != n) {
        int k = __builtin_ctz(n) - 1;
        rev.resize(n);
        for (int i = 0; i < n; ++i)
            rev[i] = rev[i >> 1] >> 1 | (i & 1) << k;
    }
    for (int i = 0; i < n; ++i)
        if (rev[i] < i)
            swap(a[i], a[rev[i]]);
    if (int(roots.size()) < n) {
        int k = __builtin_ctz(roots.size());
        roots.resize(n);
        while ((1 << k) < n) {
            std::complex<double> e = {cos(PI / (1 << k)), sin(PI / (1 << k))};
            for (int i = 1 << (k - 1); i < (1 << k); ++i) {
                roots[2 * i] = roots[i];
                roots[2 * i + 1] = roots[i] * e;
            }
            ++k;
        }
    }
    for (int k = 1; k < n; k *= 2) {
        for (int i = 0; i < n; i += 2 * k) {
            for (int j = 0; j < k; ++j) {
                auto u = a[i + j], v = a[i + j + k] * roots[k + j];
                a[i + j] = u + v;
                a[i + j + k] = u - v;
            }
        }
    }
}
void idft(std::vector<std::complex<double>> &a) {
    int n = a.size();
    reverse(a.begin() + 1, a.end());
    dft(a);
    for (int i = 0; i < n; ++i)
        a[i] /= n;
}
std::vector<ll> operator*(std::vector<ll> a, std::vector<ll> b) {
    int sz = 1, tot = a.size() + b.size() - 1;
    while (sz < tot)
        sz *= 2;
    std::vector<std::complex<double>> ca(sz), cb(sz);
    //copy(a.begin(), a.end(), ca.begin());
    //copy(b.begin(), b.end(), cb.begin());
    for(int i = 0;i < sz;i ++){
        if(i < a.size())ca[i].real(a[i]);
        if(i < b.size())ca[i].imag(b[i]);
    }
    dft(ca);
    //dft(cb);
    for (int i = 0; i < sz; ++i)
        ca[i] *= ca[i];
    idft(ca);
    a.resize(tot);
    for (int i = 0; i < tot; ++i)
        a[i] = std::floor(ca[i].imag() / 2 + 0.5);
    return a;
}
int main(){
    int n,m;
    cin >> n >> m;
    n ++,m ++;
    std::vector<ll> a(n),b(m);
    for(int i = 0;i < n;i ++)cin >> a[i];
    for(int i = 0;i < m;i ++)cin >> b[i];
    std::vector<ll> c = a * b;
    for(int i = 0;i < n + m - 1;i ++)cout << c[i] << " ";
    return 0;
}