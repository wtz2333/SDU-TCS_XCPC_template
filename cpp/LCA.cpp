struct Tree{
    int n;// 5e5
    vector<vector<int>> e;
    vector<int> dep;
    vector<array<int,21>> f;
    Tree(int n):n(n),e(n + 1),dep(n + 1),f(n + 1){}
    void clear(){
        for(int i = 0;i <= n;i ++)e[i].clear();
        dep.clear();
        f.clear();
    }
    void dfs(int x){
        for(auto v:e[x]){
            if(v == f[x][0])continue;
            f[v][0] = x;
            dep[v] = dep[x] + 1;
            dfs(v);
        }
    }
    void init(int root = 1){
        dfs(root);
        for(int j = 1;j <= 19;j ++){
            for(int i = 1;i <= n;i ++){
                f[i][j] = f[f[i][j - 1]][j - 1];
            }
        }
    }
    int lca(int u,int v){
        if(dep[u] < dep[v])swap(u,v);
        for(int i = 0;i <= 19;i ++){
            if((dep[u] - dep[v]) & (1 << i))u = f[u][i];
        }
        if(u == v)return u;
        for(int j = 19;j >= 0; j--){
            if(f[u][j] != f[v][j]){
                u = f[u][j];
                v = f[v][j];
            }
        }
        return f[u][0];
    }
    void addEdge(int u, int v) {
        e[u].push_back(v);
        e[v].push_back(u);
    }
    int kth(int x,int k){
        for(int i = 0;i <= 19;i ++){
            if(k & (1 << i))x = f[x][i];
        }
        return x;
    }
};