#include<bits/stdc++.h>

using namespace std;

typedef long long ll;

const int maxn = 3e6 + 7;
int mu[maxn],p[maxn],vis[maxn],tot;
int sum_mu[maxn];
int phi[maxn];
ll sum_phi[maxn];
map <ll,ll> mp_mu,mp_phi;
void init(const int n){
    mu[1] = 1;phi[1] = 1;
    for(int i = 2;i <= n;i ++){
        if(!vis[i]){
            p[++ tot] = i;
            mu[i] = -1;phi[i] = i - 1;
        }
        for(int j = 1;j <= tot && i * p[j] <= n;j ++){
            vis[i * p[j]] = 1;
            if(i % p[j] == 0){
            	phi[i * p[j]] = phi[i] * p[j];
                mu[i * p[j]] = 0;
            	break;
            }
            mu[i * p[j]] = -mu[i];
            phi[i * p[j]] = phi[i] * phi[p[j]];
        }
    }
    for(int i = 1;i < n;i ++)sum_mu[i] = sum_mu[i - 1] + mu[i],sum_phi[i] = sum_phi[i - 1] + phi[i];
}
ll calc_mu(ll x){
	if(x < maxn)return sum_mu[x];
	if(mp_mu[x]) return mp_mu[x];
	ll l = 2,r;ll ans = 1;
    while(l <= x){
        r = x / (x / l);
        ans -= 1ll * (r - l + 1) * calc_mu(x / l);
        l = r + 1;
    }
    return mp_mu[x] = ans;
}
ll calc_phi(ll x){
	if(x < maxn)return sum_phi[x];
	if(mp_phi[x]) return mp_phi[x];
	ll l = 2,r;ll ans = 1ll * x * (x + 1) >> 1;
    while(l <= x){
        r = x / (x / l);
        ans -= 1ll * (r - l + 1) * calc_phi(x / l);
        l = r + 1;
    }
    return mp_phi[x] = ans;
}
int main(){
	ios::sync_with_stdio(false);
	cin.tie(0);
	int T;
	cin >> T;
	init(maxn);
	while(T --){
		int n;cin >> n;
		cout << calc_phi(n) << " " << calc_mu(n) << "\n";
	}
}