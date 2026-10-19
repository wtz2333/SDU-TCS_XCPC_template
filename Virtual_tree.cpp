const int maxn = 2e5 + 7;
int dfn[maxn],mn[21][maxn];
int n,cnt;
vector<int> e[maxn],g[maxn];
int get(int x, int y) {return dfn[x] < dfn[y] ? x : y;}
void dfs(int x, int fa) {
  mn[0][dfn[x] = ++cnt] = fa;
  for(int v : e[x]) if(v != fa) dfs(v, x); 
}
int lca(int u, int v) {
  if(u == v) return u;
  if((u = dfn[u]) > (v = dfn[v])) swap(u, v);
  int d = __lg(v - u++);
  return get(mn[d][u], mn[d][v - (1 << d) + 1]);
}
void add(int x,int y){
    g[x].push_back(y);
}
int build_virtual_tree(vector<int> &h) {
  vector<int> a;
  sort(h.begin(), h.end(),[&](int &a,int &b){
      return dfn[a] < dfn[b];  
  });  // 把关键点按照 dfn 序排序
  for (int i = 0; i < h.size(); ++i) {
    a.push_back(h[i]);
    if(i + 1 != h.size())a.push_back(lca(h[i], h[i + 1]));  // 插入 lca
  }
  sort(a.begin(), a.end(), [&](int &a,int &b){
      return dfn[a] < dfn[b];  
  });  // 把所有虚树上的点按照 dfn 序排序
  a.erase(unique(a.begin(),a.end() ),a.end());
  for (int i = 0; i < a.size() - 1; ++i) {
    int lc = lca(a[i], a[i + 1]);
    add(lc, a[i + 1]);  // 连边，如有边权 就是 distance(lc,a[i+1])
  }
  return a[0];
}

int main(){
    ios::sync_with_stdio(false);
    cin.tie(0);
    int n;
    cin >> n;
    for(int i = 1;i < n;i ++){
        int a,b;
        cin >> a >> b;
        e[a].push_back(b);
        e[b].push_back(a);
    }
    dfs(1, 0);
    for(int i = 1; i <= __lg(n);i ++){
        for(int j = 1; j + (1 << i) - 1 <= n;j ++)  
        mn[i][j] = get(mn[i - 1][j], mn[i - 1][j + (1 << i - 1)]);
    }
}
