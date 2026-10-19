typedef long long ll;
const int mo = 1e9 + 7;
using i64 = long long;

template<typename T>
struct Gauss{
    int n,m;
    std::vector<std::vector<T>> a;
    Gauss(int n,int m,const vector<vector<T>> &a):n(n),m(m),a(a){}
    vector<double> gauss(){
        std::vector<double> sol;
        for(int i = 0;i < m;i ++){
            int mx = i;//最大数的行号 
            for(int j = i + 1;j < n;j ++)
                if(fabs(a[mx][i]) < fabs(a[j][i])) mx = j;
            if(fabs(a[mx][i]) < eps){
                return sol;
            }
            if(mx != i) swap(a[mx],a[i]);
            double t = a[i][i];
            for(int j = i;j < m + 1;j ++){
                a[i][j] /= t;
            }
            for(int j = i + 1;j < n;j ++){
                t = a[j][i];
                for(int k = i;k < m + 1;k ++){
                    a[j][k] -= a[i][k] * t;
                }
            }
        }
        for(int i = m;i < n;i ++){
            if(fabs(a[i][m]) > eps)return sol;
        }
        sol.resize(m);
        sol[m - 1] = a[m - 1][m];
        for(int i = n - 2;i >= 0;i --){
            sol[i] = a[i][m];
            for(int j = i + 1;j < m;j ++){
                sol[i] -= a[i][j] * sol[j];
            }
        }
        return sol;
    }
    int gauss_det(){
        assert(n == m);
        int det = 1;
        for(int i = 0;i < n;i ++){
            for(int j = i + 1;j < n;j ++){
                while(a[j][i] != 0){
                    int k = a[i][i] / a[j][i];
                    for(int h = i;h < n;h ++){
                        a[i][h] -= a[j][h] * k;
                        std::swap(a[i][h],a[j][h]);
                    }
                    det = -det;
                }
            }
        }
        for(int i = 0;i < n;i ++) det *= a[i][i];
        return det;
    }

    T gauss_det_mod(int st = 0){
        T det = T(1);
        for(int i = st;i < n;i ++){
            det *= a[i][i];
            //cerr << i << " " << a[i][i].a << " " << a[i][i].b << "\n";
            T tmp = T(1) / a[i][i]; 
            for(int j = i;j < n;j ++) a[i][j] *= tmp;
            for(int j = i + 1;j < n;j ++){
                tmp = a[j][i];
                for(int k = i;k < n;k ++) a[j][k] = a[j][k] - tmp * a[i][k];
            }
        }
        return det;
    }
};