#include<bits/stdc++.h>

using namespace std;

using ll = long long;
using pii = pair<int,int>;

const int maxn = 2e5 + 7;
int n,q;
ll w;
vector<pair<int,ll>> e[maxn];

ll c[maxn];
void add(int x,ll y) {
    for(;x <= n;x += (x & (-x))) c[x] += y;
}
void change(int l,int r,ll x) { // [l,r]
    add(l,x);add(r + 1,-x);
}
ll ask(int x) {
    ll ans = 0;
    for(;x > 0;x -= (x & (-x))) ans += c[x];
    return ans;
}

int dfn[maxn],mi[22][maxn],id[maxn],siz[maxn],cnt;
int get(int x, int y) {return dfn[x] < dfn[y] ? x : y;}
void dfs(int x, int fa) {
    mi[0][dfn[x] = ++cnt] = fa;
    id[cnt] = x;siz[x] = 1;
    for(auto [v,w] : e[x]) {
        if(v == fa) continue;
        dfs(v,x);
        change(dfn[v],dfn[v] + siz[v] - 1,w);
        siz[x] += siz[v];
    }
}
int lca(int u, int v) {
    if(u == v) return u;
    if((u = dfn[u]) > (v = dfn[v])) swap(u, v);
    int d = __lg(v - u++);
    return get(mi[d][u], mi[d][v - (1 << d) + 1]);
}
ll dis(pii a) {
    auto [u,v] = a;
    return ask(dfn[u]) + ask(dfn[v]) - 2 * ask(dfn[lca(u,v)]);
}

pii tr[maxn << 2];
pii merge(pii &a,pii &b) {
    pii p[6] = {a,b,{a.first,b.first},{a.first,b.second},{a.second,b.first},{a.second,b.second}};
    vector<ll> d(6);
    for(int i = 0;i < 6;i ++) d[i] = dis(p[i]);
    return p[max_element(d.begin(),d.end()) - d.begin()];
}
void build(int x,int l,int r) {
    if(l == r) {
        tr[x] = {id[l],id[l]};
        return ;
    }
    int mid = l + r >> 1;
    build(x << 1,l,mid);
    build(x << 1 | 1,mid + 1,r);
    tr[x] = merge(tr[x << 1],tr[x << 1 | 1]);  
}
void update(int x,int L,int R,int l,int r) {
    if(L == l && r == R) return ;
    int mid = L + R >> 1;
    if(r <= mid) update(x << 1,L,mid,l,r);
    else if(l > mid) update(x << 1 | 1,mid + 1,R,l,r);
    else update(x << 1,L,mid,l,mid),update(x << 1 | 1,mid + 1,R,mid + 1,r);
    tr[x] = merge(tr[x << 1],tr[x << 1 | 1]);  
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(0);
    cin >> n >> q >> w;
    vector<tuple<int,int,int>> edges;
    for(int i = 0;i < n - 1;i ++) {
        int x,y;ll z;
        cin >> x >> y >> z;
        e[x].push_back({y,z});
        e[y].push_back({x,z});
        edges.emplace_back(x,y,z);
    }
    dfs(1, 0);
    for(int i = 1; i <= __lg(n); i++)
      for(int j = 1; j + (1 << i) - 1 <= n; j++)
        mi[i][j] = get(mi[i - 1][j], mi[i - 1][j + (1 << i - 1)]);
    ll lastans = 0;
    build(1,1,n);
    for(int i = 0;i < q;i ++) {
        int x;ll y;
        cin >> x >> y;
        x = (x + lastans) % (n - 1);
        y = (y + lastans) % w;
        auto [son,fa,ww] = edges[x];
        if(dfn[son] < dfn[fa]) swap(son,fa);
        change(dfn[son],dfn[son] + siz[son] - 1,y - ww);
        update(1,1,n,dfn[son],dfn[son] + siz[son] - 1);
        lastans = dis(tr[1]);
        cout << lastans << "\n";
        ww = y;
        edges[x] = {son,fa,ww};
    }
}