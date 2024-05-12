# 数学

## 数论

### 扩展欧几里得（线性同余方程，斐蜀定理）

扩展欧几里得：$\gcd(a,b) = \gcd(b,a\%b)$ ，$ax + by = b{ x} + (a- \left \lfloor \frac{a}{b} \right \rfloor){y}$

斐蜀定理：$ax + by = c$ 若有解，则有 $(a,b)|c$

线性同余方程：$ax \equiv  c\pmod{b}\Rightarrow   ax + by = c$

```cpp
ll exgcd(ll a,ll b,ll &x,ll &y){
    if(b == 0){
        x = 1 ,y = 0;return a;
    }
    ll d = exgcd(b,a % b,x,y);
    ll tmp = x;
    x = y;
    y = tmp - (a / b) * y;
    return d;
}
void solve(){
    ll a,b,c;
    cin >> a >> b >> c;
    ll x0,y0;
    ll d = exgcd(a,b,x0,y0);
    if(c % d){
        cout << -1 << "\n";
        return ;
    }
    ll p = a / d,q = b / d;
    ll x = ((c / d) % q * x0 % q + q) % q;
    if(x == 0)x = q;
    ll y = (c - a * x) / b;
    if(y <= 0){
        y = ((c / d) % p * y0 % p + p) % p;
        cout << (x == 0 ? q : x) << " " << (y == 0 ? p : y) << "\n";
        return ;
    }
    ll ans_x_mn = x;
    ll ans_y_mx = y;
    y = ((c / d) % p * y0 % p + p) % p;
    if(y == 0)y = p;
    x = (c - b * y) / a;
    ll ans_x_mx = x;
    ll ans_y_mn = y;
    ll sum = min((ans_x_mx - ans_x_mn) / q,(ans_y_mx - ans_y_mn) / p);
    cout << sum + 1 << " " << ans_x_mn << " " << ans_y_mn << " " << ans_x_mx << " " << ans_y_mx << "\n";
    // 正整数解总数
}
```

### 费马小定理（逆元）

若 $p$ 为素数，$\gcd(a,p) = 1$，则 $a^{p - 1}\equiv 1\pmod{p}$

### 线性求逆元

```cpp
inv[0] = inv[1] = 1;
for(int i = 2 ;i <= n ;i++) inv[i] = (p - p/i) * inv[p % i] % p;
```

### CRT（中国剩余定理）

$$
\begin{cases}
x = b_{1}(\bmod a_{1}) \\
x = b_{2}(\bmod a_{2})\\
   \cdots \\
x = b_n (\bmod a_{n})
\end{cases}
$$

若 $ a_{1} ,a_{2},\dots,a_{n}$ 两两互质：

令 $M = \prod_{1}^{n} a_{i}, m_{i}' = \frac{M}{a_{i}},t_{i}\times m_{i}'≡1(\text{mod }a_{i})$   则有 $x= \sum_{i=1}^{n}b_{i}\times m_{i}'\times t_{i}$   （此解为唯一解）

若 $ a_{1} ,a_{2},\dots,a_{n}$​ 两两不互质：

合并两个方程组 $x = a_1p + b_1 = a_2q + b_2$

则可将方程依次两两合并为 $x≡a_{1}p+b_{1}(\bmod \text{lcm} (a_{1},a_{2}))$，其中先求解 $p$，再带入求 $x$。

```cpp
ll r1 = B[1], m1 = A[1],r2,m2;
for(int i = 1;i < n ;i ++) {
    r2 = B[i + 1],m2 = A[i + 1];
    ll a = m1,b = m2,c = r2 - r1;
    ll d = exgcd(a,b,x,y);
    if(c % d) {
        cout<<0;return 0;
    }
    ll p = a / d,q = b / d;
    x = ((x * c / d) + q) % q;
    ll mod = lcm(m2,m1);
    ll x0 = (m1 * x + r1) % mod;
    r1 = x0 < 0 ? x0 + mod : x0;
    m1 = mod;
}
cout << r1 % m1 << "\n";
```

### 卢卡斯定理

$C_n^m = C_{n\bmod p}^{m \bmod p}\cdot C_{\lfloor n/p \rfloor}^{\lfloor m/p \rfloor}$，其中 $p$ 为质数。

### 原根

满足同余式 $a^n \equiv 1(\text{mod } m)$ 的最小正整数 $n$ 存在，这个 $n$ 称作 $a$ 模 $m$ 的阶，记作 $\delta_{m}(a)$，$a,a^2,\cdots,a^{\delta_m(a)}$ 模 $m$ 两两不同余。

设 $m \in \mathbf{N}^{*}$，$g\in \mathbf{Z}$. 若 $(g,m)=1$，且 $\delta_m(g)=\varphi(m)$，则称 $g$ 为模 $m$ 的原根。

即 $g$ 满足 $\delta_m(g) = \left| \mathbf{Z}_m^* \right| = \varphi(m)$. 当 $m$ 是质数时，我们有 $g^i \bmod m,\,0 \lt i \lt m$​ 的结果互不相同。

**原根判定定理**：

设 $m \geqslant 3, (g,m)=1$，则 $g$ 是模 $m$ 的原根的充要条件是，对于 $\varphi(m)$ 的每个素因数 $p$，都有 $g^{\frac{\varphi(m)}{p}}\not\equiv 1\pmod m$​。

若一个数 $m$ 有原根，则它原根的个数为 $\varphi(\varphi(m))$​，每一个原根都形如 $g^k$ 的形式，要求满足 $\gcd(k,\varphi(n)) = 1$。

**原根存在定理**：

一个数 $m$ 存在原根当且仅当 $m=2,4,p^{\alpha},2p^{\alpha}$，其中 $p$ 为奇素数，$\alpha\in \mathbf{N}^{*}$。

### 离散对数（BSGS）

$a^{x} = b(\bmod m)$

此处只解决 $m$为质数，将 $x$ 分解为 $i\times t - p$，则有 $a^{i\times t} = b\times a^{p}(\bmod m)$​

$t = \sqrt[2]{m}$ 时均摊复杂度最小

$ 0< p < t$ 枚举 $p$ 计算出每一个 $a^{p}$ 的值存入hash表

再枚举 $i$，算出 $a^{i\times t}$​​ 的值在hash表中查找

```cpp
ll bsgs(ll a,ll b,ll p){
    map <ll,ll> hsh ;hsh.clear();
    ll t = sqrt(p) + 1,j;
    b %= p;
    for(int i = 0 ; i < t ;i++){
        ll tmp = b * qp(a,i,p) % p;
        hsh[tmp] = i;
    }
    a = qp(a,t,p);
    if(a == 0){
        if(b == 0)return 1;
        else return -1;
    }
    for(int i = 0 ;i <= t ;i++){
        ll tmp = qp(a,i,p);
        if(hsh.find(tmp) == hsh.end())j = -1;
        else  j = hsh[tmp];
        if(i * t - j >=0 && j >= 0)return i*t-j;
    }
    return -1;
}
```

### 威尔逊定理

对于素数 $p$  有 $(p-1)!\equiv -1 \pmod p$ 。

### 数论分块

```cpp
ll up(ll x,ll y){
    return (x + y - 1) / y;
}
ll calc_up(ll x){
    ll l = 1,r;ll ans = 1e18;
    while(l <= x){
        ll m = up(x,l);
        if(m == 1)r = x;else r = (x - 1) / (m - 1);
        l = r + 1;
    }
    return ans;
}
ll calc_down(ll x){
    ll l = 1,r;ll ans = 1;
    while(l <= x){
        r = x / (x / l);
        l = r + 1;
    }
    return ans;
}
```

### 积性函数

**数论函数：在所有正整数上的函数被称为算术函数（数论函数）**

**加性函数：如果数论函数 f 对于任意两个互素的正整数 $p,q$ 均有 $f(pq) = f(p) + f(q)$，称为加性函数**

**积性函数：如果数论函数 $f$ 对于任意两个互素的正整数 $p,q$ 均有 $f(pq) = f(p)f(q)$，称为积性函数**

**完全积性函数：在积性函数的基础上，$p,q$ 对于任意正整数均成立，称为完全积性函数**

若 $f(x)$ 和 $g(x)$ 均为积性函数，不难证明下列函数均为积性函数：

$$
h(x) = f(x^p)
\\
h(x) = f^p(x)
\\
h(x) = \sum_{d|x}f(d)
\\
h(x) = f(x)g(x)
$$

**常见积性函数：**

- 单位函数：$\varepsilon(n)=[n=1]$。（完全积性）

- 恒等函数：$\operatorname{id}_k(n)=n^k$，$\operatorname{id}_{1}(n)$ 通常简记作 $\operatorname{id}(n)$。（完全积性）

- 常数函数：$1(n)=1$。（完全积性）

- 除数函数：$\sigma_{k}(n)=\sum_{d\mid n}d^{k}$。$\sigma_{0}(n)$ 通常简记作 $d(n)$ 或 $\tau(n)$，$\sigma_{1}(n)$ 通常简记作 $\sigma(n)$。

- 欧拉函数：$\varphi(n)=\sum_{i=1}^n[\gcd(i,n)=1]$

- 莫比乌斯函数：$\mu(n)=\begin{cases}1&n=1\\0&\exists d>1,d^{2}\mid n\\(-1)^{\omega(n)}&\text{otherwise}\end{cases}$，其中 $\omega(n)$ 表示 $n$ 的本质不同质因子个数，它是一个加性函数。

### 线性筛

**一般情况下可通过线性筛快速筛出积性函数**

```cpp
void init(const int n){
    mu[1] = 1;phi[1] = 1;
    for(int i = 2;i < n;i ++){
        if(!vis[i]){
            p[++ tot] = i;
            mu[i] = -1;phi[i] = i - 1;
        }
        for(int j = 1;j <= tot && i * p[j] < n;j ++){
            vis[i * p[j]] = 1;
            if(i % p[j] == 0){
                phi[i * p[j]] = phi[i] * p[j];
                mu[i * p[j]] = 0;
                break;
            }
            mu[i * p[j]] = -mu[i];
            phi[i * p[j]] = phi[i] * phi[p[j]];
        }
    }
}
```

### 欧拉函数

欧拉函数（Euler's totient function），即 $\varphi(n)$，表示的是小于等于 $n$ 和 $n$ 互质的数的个数。$\varphi(n) = \sum_{i = 1}^n [\gcd(i,n)=1]$

由唯一分解定理，设 $n = \prod_{i=1}^{s}p_i^{k_i}$，其中 $p_i$ 是质数，有 $\varphi(n) = n \times \prod_{i = 1}^s{\dfrac{p_i - 1}{p_i}}$。

如果 $(a,b) = 1$ , $\varphi (a*b) = \varphi (a) * \varphi (b)$
如果 $a$ 或 $b$为质数 $\varphi (a*b) = \varphi (a) * \varphi (b)$
如果$(a,b) \neq 1$, $\varphi (a*b) = \varphi (a) * b$

### 欧拉定理及扩展

如果 $(a,m),a^{\varphi (m)} \equiv 1(\text{mod } m)$
当 $b \geq \varphi(p)$ 时 $a^b \equiv a^{b \text{ mod } \varphi(p) + φ(p)} (\text{mod } p)$
当 $b < \varphi(p)$ 时 $a^b \equiv a^{b} (\text{mod } p)$

### 狄利克雷卷积

对于两个数论函数 $f(x)$ 和 $g(x)$，则它们的狄利克雷卷积得到的结果 $h(x)$ 定义为：
$h(x)=\sum_{d\mid x}{f(d)g\left(\dfrac xd \right)}=\sum_{ab=x}{f(a)g(b)}$
上式可以简记为：$h=f*g $

狄利克雷卷积满足交换律，结合律，分配律

单位函数 $\varepsilon$ 是 Dirichlet 卷积运算中的单位元，即对于任何数论函数 $f$，都有 $f*\varepsilon=f$

对于任何一个满足 $f(x)\ne 0$ 的数论函数，如果有另一个数论函数 $g(x)$ 满足 $f*g=\varepsilon$，则称 $g(x)$ 是 $f(x)$ 的逆元。由 **等式的性质** 可知，逆元是唯一的

常见数论卷积

$\phi * 1 = id$

$\mu * 1 = \varepsilon$

$\mu * id = \phi$

**两个积性函数的 Dirichlet 卷积也是积性函数**

**积性函数的逆元也是积性函数，且 $f(1) \neq 0$**

### 莫比乌斯反演

莫比乌斯函数

$$
\mu(n) =
\begin{cases}
1&n = 1\\
0&n\text{含有平方因子}\\
(-1)^k&k\text{为}n\text{的本质不同质因子}
\end{cases}
$$

$\mu * 1 = \varepsilon$

形式一：

$f(n) = \sum_{d|n}g(d)\Rightarrow g(n) = \sum_{d|n}\mu(d)f(\frac{n}{d})$

形式二：

$f(n) = \sum_{n|d}g(d)\Rightarrow g(n) = \sum_{n|d}\mu(\frac{d}{n})f(d)$

上两式的证明可通过等价替换和交换求和号来推导

**此外我们也可以通过狄利克雷卷积来理解**

$f = g * 1,g = f * \mu$

实际应用中我们常使用 $\varepsilon(\gcd) = \sum_{d|\gcd}\mu(d)$，即 $\mu * 1 = \varepsilon$

同时也是形式1中$f = \varepsilon$ 的情况

### 欧拉反演

$\varphi *1 = id$

展开形式同莫比乌斯反演，本质上是莫比乌斯反演的进一步推导，卷积式也可用 $\mu * 1 = \varepsilon$ 推出

实际应用中常使用 $\gcd = \sum_{d|\gcd}\varphi(d)$

### 杜教筛

对于数论函数 $f$，要计算 $S(n) = \sum_{i = 1}^nf(i)$

找到一个数论函数 $g$，有 $\sum_{i = 1}^n(f*g)(i) = \sum_{i = 1}^n g(i)S(\lfloor\frac{n}{i}\rfloor)$

得到 $g(1)S(n) = S(n) = \sum_{i = 1}^n(f*g)(i) - \sum_{i = 2}^n g(i)S(\lfloor \frac{n}{i}\rfloor )$

```cpp
const int maxn = 3e6 + 7;
int mu[maxn],p[maxn],vis[maxn],tot;
int sum_mu[maxn];
int phi[maxn];
ll sum_phi[maxn];
map <ll,ll> mp_mu,mp_phi;
void init(const int n){
    mu[1] = 1;phi[1] = 1;
    for(int i = 2;i <= n;i ++){
        if(!vis[i]){
            p[++ tot] = i;
            mu[i] = -1;phi[i] = i - 1;
        }
        for(int j = 1;j <= tot && i * p[j] <= n;j ++){
            vis[i * p[j]] = 1;
            if(i % p[j] == 0){
                phi[i * p[j]] = phi[i] * p[j];
                mu[i * p[j]] = 0;
                break;
            }
            mu[i * p[j]] = -mu[i];
            phi[i * p[j]] = phi[i] * phi[p[j]];
        }
    }
    for(int i = 1;i < n;i ++)sum_mu[i] = sum_mu[i - 1] + mu[i],sum_phi[i] = sum_phi[i - 1] + phi[i];
}
ll calc_mu(ll x){
    if(x < maxn)return sum_mu[x];
    if(mp_mu[x]) return mp_mu[x];
    ll l = 2,r;ll ans = 1;
    while(l <= x){
        r = x / (x / l);
        ans -= 1ll * (r - l + 1) * calc_mu(x / l);
        l = r + 1;
    }
    return mp_mu[x] = ans;
}
ll calc_phi(ll x){
    if(x < maxn)return sum_phi[x];
    if(mp_phi[x]) return mp_phi[x];
    ll l = 2,r;ll ans = 1ll * x * (x + 1) >> 1;
    while(l <= x){
        r = x / (x / l);
        ans -= 1ll * (r - l + 1) * calc_phi(x / l);
        l = r + 1;
    }
    return mp_phi[x] = ans;
}
```

由于符合数论反演实际意义的函数不多所以大部分数论反演题目基本上都是对上述两种反演的卷积式的应用，变形之后进行求和号交换，换元等数学手段处理等到可以快速求解的算式

### 素数测试与因式分解（Miller-Rabin & Pollard-Rho）

```cpp
i64 mul(i64 a, i64 b, i64 m) {
    return static_cast<__int128>(a) * b % m;
}
i64 power(i64 a, i64 b, i64 m) {
    i64 res = 1 % m;
    for (; b; b >>= 1, a = mul(a, a, m))
        if (b & 1)
            res = mul(res, a, m);
    return res;
}
bool isprime(i64 n) {
    if (n < 2)
        return false;
    static constexpr int A[] = {2, 3, 5, 7, 11, 13, 17, 19, 23};
    int s = __builtin_ctzll(n - 1);
    i64 d = (n - 1) >> s;
    for (auto a : A) {
        if (a == n)
            return true;
        i64 x = power(a, d, n);
        if (x == 1 || x == n - 1)
            continue;
        bool ok = false;
        for (int i = 0; i < s - 1; ++i) {
            x = mul(x, x, n);
            if (x == n - 1) {
                ok = true;
                break;
            }
        }
        if (!ok)
            return false;
    }
    return true;
}
std::vector<i64> factorize(i64 n) {
    std::vector<i64> p;
    std::function<void(i64)> f = [&](i64 n) {
        if (n <= 10000) {
            for (int i = 2; i * i <= n; ++i)
                for (; n % i == 0; n /= i)
                    p.push_back(i);
            if (n > 1)
                p.push_back(n);
            return;
        }
        if (isprime(n)) {
            p.push_back(n);
            return;
        }
        auto g = [&](i64 x) {
            return (mul(x, x, n) + 1) % n;
        };
        i64 x0 = 2;
        while (true) {
            i64 x = x0;
            i64 y = x0;
            i64 d = 1;
            i64 power = 1, lam = 0;
            i64 v = 1;
            while (d == 1) {
                y = g(y);
                ++lam;
                v = mul(v, std::abs(x - y), n);
                if (lam % 127 == 0) {
                    d = std::gcd(v, n);
                    v = 1;
                }
                if (power == lam) {
                    x = y;
                    power *= 2;
                    lam = 0;
                    d = std::gcd(v, n);
                    v = 1;
                }
            }
            if (d != n) {
                f(d);
                f(n / d);
                return;
            }
            ++x0;
        }
    };
    f(n);
    std::sort(p.begin(), p.end());
    return p;
}
```

### 公式

#### 一些数论公式

- 当 $x\geq\phi(p)$ 时有 $a^x\equiv a^{x \; mod \; \phi(p) + \phi(p)}\pmod p$
- $\mu^2(n)=\sum_{d^2|n} \mu(d)$
- $\sum_{d|n} \varphi(d)=n$
- $\sum_{d|n} 2^{\omega(d)}=\sigma_0(n^2)$，其中 $\omega$ 是不同素因子个数
- $\sum_{d|n} \mu^2(d)=2^{\omega(d)}$

#### 一些数论函数求和的例子

+ $\sum_{i=1}^n i[gcd(i, n)=1] = \frac {n \varphi(n) + [n=1]}{2}$
+ $\sum_{i=1}^n \sum_{j=1}^m [gcd(i,j)=x]=\sum_d \mu(d) \lfloor \frac n {dx} \rfloor  \lfloor \frac m {dx} \rfloor$
+ $\sum_{i=1}^n \sum_{j=1}^m gcd(i, j) = \sum_{i=1}^n \sum_{j=1}^m \sum_{d|gcd(i,j)} \varphi(d) = \sum_{d} \varphi(d) \lfloor \frac nd \rfloor \lfloor \frac md \rfloor$
+ $S(n)=\sum_{i=1}^n \mu(i)=1-\sum_{i=1}^n \sum_{d|i,d < i}\mu(d) \overset{t=\frac id}{=} 1-\sum_{t=2}^nS(\lfloor \frac nt \rfloor)$
  + 利用 $[n=1] = \sum_{d|n} \mu(d)$
+ $S(n)=\sum_{i=1}^n \varphi(i)=\sum_{i=1}^n i-\sum_{i=1}^n \sum_{d|i,d<i} \varphi(i)\overset{t=\frac id}{=} \frac {i(i+1)}{2} - \sum_{t=2}^n S(\frac n t)$
  + 利用 $n = \sum_{d|n} \varphi(d)$
+ $\sum_{i=1}^n \mu^2(i) = \sum_{i=1}^n \sum_{d^2|n} \mu(d)=\sum_{d=1}^{\lfloor \sqrt n \rfloor}\mu(d) \lfloor \frac n {d^2} \rfloor$
+ $\sum_{i=1}^n \sum_{j=1}^n gcd^2(i, j)= \sum_{d} d^2 \sum_{t} \mu(t) \lfloor \frac n{dt} \rfloor ^2 \\
  \overset{x=dt}{=} \sum_{x} \lfloor \frac nx \rfloor ^ 2 \sum_{d|x} d^2 \mu(\frac xd)$
+ $\sum_{i=1}^n \varphi(i)=\frac 12 \sum_{i=1}^n \sum_{j=1}^n [i \perp j] - 1=\frac 12 \sum_{i=1}^n \mu(i) \cdot\lfloor \frac n i \rfloor ^2-1$

#### 斐波那契数列性质

- $F_{a+b}=F_{a-1} \cdot F_b+F_a \cdot F_{b+1}$
- $F_1+F_3+\dots +F_{2n-1} = F_{2n},F_2 + F_4 + \dots + F_{2n} = F_{2n + 1} - 1$
- $\sum_{i=1}^n F_i = F_{n+2} - 1$
- $\sum_{i=1}^n F_i^2 = F_n \cdot F_{n+1}$
- $F_n^2=(-1)^{n-1} + F_{n-1} \cdot F_{n+1}$
- $gcd(F_a, F_b)=F_{gcd(a, b)}$
- 模 $n$ 周期（皮萨诺周期）
  - $\pi(p^k) = p^{k-1} \pi(p)$
  - $\pi(nm) = lcm(\pi(n), \pi(m)), \forall n \perp m$
  - $\pi(2)=3, \pi(5)=20$
  - $\forall p \equiv \pm 1\pmod {10}, \pi(p)|p-1$
  - $\forall p \equiv \pm 2\pmod {5}, \pi(p)|2p+2$
