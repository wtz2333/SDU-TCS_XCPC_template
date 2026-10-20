const int mo = 998244353;
int qp(int a,int b){
    int ans = 1,base = a % mo;
    while(b != 0){
        if(b & 1)ans = 1ll * ans * base % mo;
        base = 1ll * base * base % mo;
        b >>= 1;
    }
    return ans;
}
vector<int> lagrange(const vector<int> &x,const vector<int> y){
    assert(x.size() == y.size());
    int n = x.size();
    std::vector<int> a(n);
    for(int i = 0;i < n;i ++){
        int A = 1;
        for(int j = 0;j < n;j ++){
            if(i == j)continue;
            assert(x[i] - x[j]);
            A = 1ll * A * (x[i] - x[j] + mo) % mo;
        }
        a[i] = 1ll * y[i] * qp(A,mo - 2) % mo;
    }
    std::vector<int> b(n + 1),c(n),f(n);
    b[0] = 1;
    for(int i = 0;i < n;i ++){
        for(int j = i + 1;j >= 1;j --){
            b[j] = (1ll * b[j] * (mo - x[i]) % mo + b[j - 1]) % mo;
        }
        b[0] = 1ll * b[0] * (mo - x[i]) % mo;
    }
    for(int i = 0;i < n;i ++){
        int inv = qp(mo - x[i],mo - 2);
        if(!inv){
            for(int j = 0;j < n;j ++)c[j] = b[j + 1];
        }else {
            c[0] = 1ll * b[0] * inv % mo;
            for(int j = 1;j < n;j ++){
                c[j] = 1ll * (b[j] - c[j - 1] + mo) * inv % mo;
            }
        }
        for(int j = 0 ;j < n;j ++){
            f[j] = (f[j] + 1ll * a[i] * c[j] % mo) % mo;
        }
    }
    return f;
}