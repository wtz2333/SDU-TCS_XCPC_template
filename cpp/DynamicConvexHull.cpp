struct Item {
    P p;
    mutable P vec;
    int q = 0;
};

bool operator<(const Item &a, const Item &b) {
    if (!b.q) {
        return a.p.x < b.p.x;
    }
    return dot(a.vec, b.p) > 0;
}

struct Hull {
    std::set<Item> s;
    i128 dx = 0;
    i128 dy = 0;
};

void print(const Hull &h) {
    for (auto it : h.s) {
        std::cerr << "(" << i64(it.p.x + h.dx) << ", " << i64(it.p.y + h.dy) << ") ";
    }
    std::cerr << "\n";
}

constexpr i64 inf = 2E18;

void insert(Hull &h, P p) {
    p.x -= h.dx;
    p.y -= h.dy;
    h.s.insert({p});
    auto it = h.s.lower_bound({p});
    if (it != h.s.end() && it->p.x == p.x) {
        if (it->p.y > p.y) {
            return;
        }
        it = h.s.erase(it);
    }
    if (it != h.s.begin() && it != h.s.end()
        && cross(p - std::prev(it)->p, it->p - p) >= 0) {
        return;
    }
    it = h.s.insert({p}).first;
    auto r = std::next(it);
    if (r != h.s.end()) {
        while (cross(r->p - p, r->vec) >= 0) {
            r = h.s.erase(r);
        }
        it->vec = r->p - p;
    } else {
        it->vec = P(0, -inf);
    }
    
    if (it != h.s.begin()) {
        auto l = std::prev(it);
        while (l != h.s.begin()) {
            auto a = std::prev(l);
            if (cross(a->vec, p - l->p) < 0) {
                break;
            }
            h.s.erase(l);
            l = a;
        }
        l->vec = p - l->p;
    }
}

i64 query(const Hull &h, i64 x) {
    if (h.s.empty()) {
        return 0LL;
    }
    auto it = h.s.lower_bound({P(x, 1), P{}, 1});
    assert(it != h.s.end());
    auto p = it->p;
    p.x += h.dx;
    p.y += h.dy;
    return p.x * x + p.y;
}