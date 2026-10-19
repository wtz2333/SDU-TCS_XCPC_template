constexpr long long INF = 1'000'000'000'000'000'000;
constexpr int C = 100'000;
struct Line {
    int k;
    long long b;
    Line(int k, long long b) : k(k), b(b) {}
};
long long f(const Line &line, int x) {
    return 1LL * line.k * x + line.b;
}
struct Node {
    Node *lc, *rc;
    Line line;
    Node(const Line &line) : lc(nullptr), rc(nullptr), line(line) {}
};
void modify(Node *&p, int l, int r, Line line) {
    if (p == nullptr) {
        p = new Node(line);
        return;
    }
    int m = (l + r) / 2;
    bool le = f(p -> line, l) < f(line, l);
    bool mi = f(p -> line, m) < f(line, m);
    if (!mi)
        std::swap(p -> line, line);
    if (r - l == 1)
        return;
    if (le != mi) {
        modify(p -> lc, l, m, line);
    } else {
        modify(p -> rc, m, r, line);
    }
}
Node *merge(Node *p, Node *q, int l, int r) {
    if (p == nullptr)
        return q;
    if (q == nullptr)
        return p;
    int m = (l + r) / 2;
    p -> lc = merge(p -> lc, q -> lc, l, m);
    p -> rc = merge(p -> rc, q -> rc, m, r);
    modify(p, l, r, q -> line);
    return p;
}
long long query(Node *p, int l, int r, int x) {
    if (p == nullptr)
        return INF;
    long long ans = f(p -> line, x);
    if (r - l == 1)
        return ans;
    int m = (l + r) / 2;
    if (x < m) {
        return std::min(ans, query(p -> lc, l, m, x));
    } else {
        return std::min(ans, query(p -> rc, m, r, x));
    }
}