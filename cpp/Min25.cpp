#include<bits/stdc++.h>

using namespace std;
using ll = long long;
const int mo = 1e9 + 7;
std::vector<int> minp, primes;
void sieve(int n) {
    minp.assign(n + 1, 0);
    primes.clear();
    for (int i = 2; i <= n; i++) {
        if (minp[i] == 0) {
            minp[i] = i;
            primes.push_back(i);
        }
        for (auto p : primes) {
            if (i * p > n) break;
            minp[i * p] = p;
            if (p == minp[i]) {
                break;
            }
        }
    }
}
int qp(int a,int b) {
    int ans = 1,base = a;
    while(b != 0) {
        if(b & 1) ans = 1ll * ans * base % mo;
        base = 1ll * base * base % mo;
        b >>= 1;
    }
    return ans;
}
int main() {
    ios::sync_with_stdio(false);
    cin.tie(0);
    ll n;
    cin >> n;
    ll s = sqrt(n);
    sieve(s);
    vector<ll> v;
    for(ll l = 1,r;l <= n;l = r + 1) {
        v.push_back(n / l);
        r = n / (n / l);
    }
    int inv2 = qp(2,mo - 2);
    int inv6 = qp(6,mo - 2);
    auto calc1 = [&](ll x) ->ll {
        x %= mo;
        return (x * (x + 1) % mo) * inv2 % mo;
    };
    auto calc2 = [&](ll x) ->ll {
        x %= mo;
        return (x * (x + 1) % mo * (2 * x + 1) % mo) * inv6 % mo;
    };
    vector<ll> f1(v.size() + 1),f2(v.size() + 1);
    auto getid = [&](ll x) -> ll {
        if(x <= s) return v.size() - x;
        else return n / x - 1;
    };
    for(int i = 0;i < v.size();i ++) {
        ll x = v[i];
        x %= mo;
        f1[i] = (calc1(x) + mo - 1) % mo;// \sum_2^n i^k = \sum_1^n i^k - 1
        f2[i] = (calc2(x) + mo - 1) % mo;
    }
    ll ps1 = 0,ps2 = 0;
    for(auto p : primes) {
        for(int i = 0;i < v.size();i ++) {
            if(1ll * p * p > v[i]) break;
            f1[i] -= 1ll * p * (f1[getid(v[i] / p)] - ps1 + mo) % mo;
            f1[i] += mo;f1[i] %= mo;
            f2[i] -= 1ll * p * p % mo * (f2[getid(v[i] / p)] - ps2 + mo) % mo;
            f2[i] += mo;f2[i] %= mo;
        }
        ps1 += p;ps1 %= mo;
        ps2 += 1ll * p * p % mo;ps2 %= mo;
    }
    auto F = [](ll x) -> ll { // targeted function
        x %= mo;
        return (x * x % mo - x + mo) % mo;
    };
    auto S = [&](auto self,ll x,int y) -> ll {
        int p = y == 0 ? 0 : primes[y - 1];
        if(p >= x) return 0ll;
        ll res = (f2[getid(x)] - f2[getid(p)] - (f1[getid(x)] - f1[getid(p)]) + mo + mo) % mo;
        for(int i = y;i < primes.size() && primes[i] <= x / primes[i];i ++) {
            ll w = primes[i];
            for(int j = 1;w <= x;j ++,w = w * primes[i]) {
                res = (res + F(w) * (self(self,x / w,i + 1) % mo + (j != 1)) % mo) % mo;
            }
        }
        return res;
    };
    cout << (S(S,n,0) + 1) % mo << "\n";
}