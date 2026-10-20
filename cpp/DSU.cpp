struct DSU {
    std::vector<int> f, siz;
    DSU(int n) : f(n), siz(n, 1) { std::iota(f.begin(), f.end(), 0); }
    int leader(int x) {
        while (x != f[x]) x = f[x] = f[f[x]];
        return x;
    }
    bool same(int x, int y) { return leader(x) == leader(y); }
    bool merge(int x, int y) {
        x = leader(x);
        y = leader(y);
        if (x == y) return false;
        siz[x] += siz[y];
        f[y] = x;
        return true;
    }
    int size(int x) { return siz[leader(x)]; }
};

struct DSU {
    std::vector<std::pair<int &, int>> his;
    
    int n,c,o;
    std::vector<int> f, g;
    
    DSU(int n_) : n(n_), f(n + 1, -1), g(n + 1), c(n), o(0) {}
    
    std::pair<int, int> find(int x) {
        if (f[x] < 0) {
            return {x, 0};
        }
        auto [u, v] = find(f[x]);
        return {u, v ^ g[x]};
    }
    
    void set(int &a, int b) {
        his.emplace_back(a, a);
        a = b;
    }
    
    void merge(int a, int b) {
        auto [u, xa] = find(a);
        auto [v, xb] = find(b);
        int w = xa ^ xb ^ 1;
        if (u == v) {
            if(w) set(o,o + 1);
            return;
        }
        set(c,c - 1);
        if (f[u] > f[v]) {
            std::swap(u, v);
        }
        set(f[u], f[u] + f[v]); // this f = siz
        set(f[v], u);
        set(g[v], w);
    }
    
    int timeStamp() {
        return his.size();
    }
    
    void rollback(int t) {
        while (his.size() > t) {
            auto [x, y] = his.back();
            x = y;
            his.pop_back();
        }
    }
};