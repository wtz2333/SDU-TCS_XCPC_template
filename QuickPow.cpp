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