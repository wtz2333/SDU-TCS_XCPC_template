struct Persistent_SegmentTree {
    int sum[(maxn << 5) + 10], root[maxn + 10], lch[(maxn << 5) + 10],
        rch[(maxn << 5) + 10], cnt[(maxn << 5) + 10];
    int tot = 0;
    void init() {tot = 0;}
    int a[maxn + 10];
    void update(int &rot, int pr, int L, int R, int k) { // 插入操作
        rot = ++tot;
        lch[rot] = lch[pr];
        rch[rot] = rch[pr];
        sum[rot] = sum[pr] + k;
        cnt[rot] = cnt[pr] + 1;
        if (L == R)
            return;
        int mid = (L + R) >> 1;
        if (k <= mid)
            update(lch[rot], lch[pr],  L, mid,k);
        else
            update(rch[rot], rch[pr],  mid + 1, R,k);
    }

    int getsum(int s, int t, int L, int R, int l, int r) { // s,t为root[s],root[t]的根节点 中所有大小在[l,r]之间数字的和
        if (l <= L && R <= r)
            return sum[t] - sum[s];
        int res = 0;
        int mid = (L + R) >> 1;
        if (l <= mid)
            res += getsum(lch[s], lch[t], L, mid, l, r);
        if (r > mid)
            res += getsum(rch[s], rch[t], mid + 1, R, l, r);
        return res;
    }
    int get_Kth_min_Sum(int s,int t,int l,int r,int k,int &nowsum,int &ans) { //return 第k小的值
        int ss = nowsum + cnt[t] - cnt[s];
        if (ss < k) {
            nowsum = ss;
            ans += sum[t] - sum[s];
            return -1;
        }
        if (l == r){
            ans += (k - nowsum) * l;
            return l;
        }
        int mid = (l + r) / 2;
        int pos = get_Kth_min_Sum(lch[s],lch[t],l,mid,k,nowsum,ans);
        if (pos != -1)return pos;
        return get_Kth_min_Sum(rch[s],rch[t],mid + 1,r,k,nowsum,ans);
    }
    int get_Kth_max_Sum(int s,int t,int l,int r,int k,int &nowsum,int &ans) { //return 第k大的值
        int ss = nowsum + cnt[t] - cnt[s];
        if (ss < k) {
            nowsum = ss;
            ans += sum[t] - sum[s];
            return -1;
        }
        if (l == r){
            ans += (k - nowsum) * l;
            return l;
        }
        int mid = (l + r) / 2;
        int pos = get_Kth_max_Sum(rch[s],rch[t],mid + 1,r,k,nowsum,ans);
        if (pos != -1)return pos;
        return get_Kth_max_Sum(lch[s],lch[t],l,mid,k,nowsum,ans);
    }
    int get_upper(int s,int t,int l,int r,int x) { //第一个大于等于x的数
        if (r < x)return -1;
        if (x <= l) {
            int ss = cnt[t] - cnt[s];
            if (ss == 0) {
                return -1;
            }
            if (l == r)return l;
        }
        int mid = (l + r) / 2;
        int pos = get_upper(lch[s],lch[t],l,mid,x);
        if (pos != -1)return pos;
        return get_upper(rch[s],rch[t],mid + 1,r,x);
    }
    int get_lower(int s,int t,int l,int r,int x) { //第一个小于等于x的数
        if(l > x)return -1;
        if(x >= r){
            int ss = cnt[t] - cnt[s];
            if(ss == 0){
                return -1;
            }
            if(l == r)return r;
        }
        int mid = (l + r) / 2;
        int pos = get_lower(rch[s],rch[t],mid + 1,r,x);
        if(pos != -1)return pos;
        return get_lower(lch[s],lch[t],l,mid,x);
    }
};