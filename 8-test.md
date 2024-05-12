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


设 $m \in \mathbf{N}^{*}$，$g\in \mathbf{Z}$. 若 $(g,m)=1$，且 $\delta_m(g)=\varphi(m)$，则称 $g$ 为模 $m$ 的原根。

即 $g$ 满足 $\delta_m(g) = \left| \mathbf{Z}_m^{*} \right| = \varphi(m)$. 当 $m$ 是质数时，我们有 $g^i \bmod m,\,0 < i < m$​ 的结果互不相同。

**原根判定定理**：

设 $m \ge 3, (g,m)=1$，则 $g$ 是模 $m$ 的原根的充要条件是，对于 $\varphi(m)$ 的每个素因数 $p$，都有 $g^{\frac{\varphi(m)}{p}}\not\equiv 1\pmod m$​。

若一个数 $m$ 有原根，则它原根的个数为 $\varphi(\varphi(m))$​，每一个原根都形如 $g^k$ 的形式，要求满足 $\gcd(k,\varphi(n)) = 1$。
