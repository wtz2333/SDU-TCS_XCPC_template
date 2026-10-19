struct Hash {
	string s;
	using ull = unsigned long long;
	ull P1 = 998255347;
	ull P2 = 1018253347;
	ull base = 131;
	vector<ull> hs1,hs2;
	vector<ull> ps1,ps2;
	Hash(string s): s(s) {
		int n = s.size();
		hs1.resize(n);
		hs2.resize(n);
		ps1.resize(n);
		ps2.resize(n);
		ps1[0] = ps2[0] = 1;
		hs1[0] = hs2[0] = (s[0] - 'a' + 1);
		for(int i = 1;i < n;i ++) {
			hs1[i] = hs1[i - 1] * base % P1 + (s[i] - 'a' + 1);
			hs2[i] = hs2[i - 1] * base % P2 + (s[i] - 'a' + 1);
			ps1[i] = (ps1[i - 1] * base) % P1;
			ps2[i] = (ps2[i - 1] * base) % P2;
		}
	}
	pair<ull,ull> query(int l,int r) {
		ull res1 = (hs1[r] - (l == 0 ? 0 : hs1[l - 1]) * ps1[r - l + 1] % P1 + P1) % P1;
		ull res2 = (hs2[r] - (l == 0 ? 0 : hs2[l - 1]) * ps2[r - l + 1] % P2 + P2) % P2;
		return {res1,res2};
	} // [l,r]
};