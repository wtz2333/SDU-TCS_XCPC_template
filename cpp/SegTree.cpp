struct SegTree {
	int l, r;
	SegTree *ls, *rs;
	ll sum;
	ll plus;
    SegTree (const int L, const int R) : l(L), r(R) {
		plus = 0;
		if (L == R) {
			/*Initial*/
			sum = 0;
			ls = rs = nullptr;
		} else {
			int M = (L + R) >> 1;
			ls = new SegTree (L, M);
			rs = new SegTree (M + 1, R);
			pushup();
		}
	}
	void pushup() {
		sum = ls->sum + rs->sum;
		//  std::cerr << "AAA" << l << ' ' << r << ' ' << sum;
	}
	void make_tag(long long w) {
		sum += (r - l + 1) * w;
		plus += w;
	}
	void pushdown() {
		if (plus == 0) return;
		ls->make_tag(plus);
		rs->make_tag(plus);
		plus = 0;
	}
	void upd(const int L, const int R, const int w) {
		if ((L > r) || (l > R)) return;
		if ((L <= l) && (r <= R)) {
			make_tag(w);
		} else {
			pushdown();
			ls->upd(L, R, w);
			rs->upd(L, R, w);
			pushup();
		}
	}
	void upd(const int x,const int w) {
        if((x > r) || (l > x)) return ;
        if(l == x && r == x) {
            sum = w;
        } else {
            ls->upd(x,w);
            rs->upd(x,w);
            pushup();
        }
    }
    int qry(const int L,const int R) {
        if ((L > r) || (l > R)) return 0;
        if ((L <= l) && (r <= R)) {
            return sum;
        } else {
            pushdown();
            return ls->qry(L, R) + rs->qry(L, R);
        }
    }
};