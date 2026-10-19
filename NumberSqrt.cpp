ll up(ll x,ll y){
    return (x + y - 1) / y;
}
ll calc_up(ll x){   
    ll l = 1,r;ll ans = 1e18;
    while(l <= x){
        ll m = up(x,l);
        if(m == 1)r = x;else r = (x - 1) / (m - 1);
        l = r + 1;
    }
    return ans;
}
ll calc_down(ll x){
    ll l = 1,r;ll ans = 1;
    while(l <= x){
        r = x / (x / l);
        l = r + 1;
    }
    return ans;
}