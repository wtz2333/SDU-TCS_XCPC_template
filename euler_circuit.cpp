struct Graph{
    int n;
    struct Edge {
        int to;
        Edge(int to) : to(to) {}
    };
    vector<Edge> e;
    vector<vector<int>> g;
    vector<int> cur,s;
    vector <int> du;
    Graph(int n) : n(n) , g(n){
        cur.assign(n,0);
        du.assign(n,0);
    }
    void init(int n) {
        for (int i = 0; i < n; i++) g[i].clear();
        e.clear();s.clear();
    }
    void add(int x,int y){
        g[x].push_back(e.size());
        e.push_back(y);
        du[y] ++;
    }
    void dfs(int u) {
        //cout << u << "\n";
        for(int i = cur[u];i < g[u].size();i = cur[u]){
            int j = g[u][i];
            int v = e[j].to;
            cur[u] = i + 1;
            dfs(v);
        }
        s.push_back(u);
    }
    void print(){
        for(int i = s.size() - 1;i >= 0;i --)cout << s[i] << (i == 0 ? "\n" : " "); 
    }
};