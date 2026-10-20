struct LinearBasis {
    static const int n = 60;
    bool O; ll a[n + 1]; int cnt;

    void clear() { O = 0; cnt = 0; for (int i = 0; i <= n; ++i) a[i] = 0; } 
    void insert(ll v) {
        for (int i = n; ~i; --i) {
            if (!(v >> i & 1)) continue;
            if (!a[i]) { a[i] = v; break; }
            v ^= a[i]; 
        } O |= !v; 
    }
    void work() {
        for (int i = 0; i <= n; cnt += a[i++] > 0)
            for (int j = 0; j < i; ++j)
                if (a[i] >> j & 1) a[i] ^= a[j];
    }
    ll get_max(ll v = 0) {
        for (int i = n; ~i; --i) v = max(v, v ^ a[i]); 
        return v;
    }
    ll get_min() {
        if (O) return 0;
        for (int i = 0; i <= n; ++i)
            if (a[i]) return a[i];
    }
    ll get_min(ll v) {
        for (int i = n; ~i; --i) v = min(v, v ^ a[i]);
        return v;
    } 
    ll get_kth(ll k) {
        ll ans = 0; k -= O; 
        for (int i = 0; i <= n; ++i) {
            if (!a[i]) continue;
            if (k & 1) ans ^= a[i];
            k >>= 1; 
        }
        if (k > 0) return -1; 
        else return ans; 
    }
    bool check(ll v) {
        for (int i = n; ~i; --i) {
            if (!(v >> i & 1)) continue;
            if (!a[i]) return 0;
            v ^= a[i]; 
        } return 1;
    } 
} S;