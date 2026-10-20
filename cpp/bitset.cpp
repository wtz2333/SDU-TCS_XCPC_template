#include<vector>
using ull = unsigned long long;
struct Bit {
    ull mi[65];
    // ull bit[15626];
    std::vector<ull> bit; int len;
    Bit() {
        len = 10;
        bit.resize(len);
        for(int i = 0;i <= 63;i ++) mi[i] = (1ull << i);
    }
    Bit(int len) : len(len){
        bit.resize(len);
        for(int i = 0;i <= 63;i ++) mi[i] = (1ull << i);
    }
    void reset() {bit.assign(len,0);}
    void set1(int x) { bit[x>>6] |= mi[x&63];}
    void set0(int x) { bit[x>>6] &= ~mi[x&63];}
    void flip(int x) { bit[x>>6] ^= mi[x&63];}
    bool operator [](int x) { 
        return (bit[x>>6] >> (x&63)) & 1;
    }   
    int count() {
        int s = 0;
        for(int i = 0;i < len;i ++)s += __builtin_popcountll(bit[i]);
        return s;
    } 
    Bit operator ~ (void) const {
        Bit res;
        for (int i = 0;i < len; i++) res.bit[i] = ~bit[i];
        return res;
    }
    
    Bit operator & (const Bit &b) const {
        Bit res;
        for (int i = 0;i < len; i++) res.bit[i] = bit[i] & b.bit[i];
        return res;
    }
    
    Bit operator | (const Bit &b) const {
        Bit res;
        for (int i = 0;i < len; i++) res.bit[i] = bit[i] | b.bit[i];
        return res;
    }   
    
    Bit operator ^ (const Bit &b) const {
        Bit res;
        for (int i = 0;i < len; i++) res.bit[i] = bit[i] ^ b.bit[i];
        return res;
    }
    
    void operator &= (const Bit &b)  {
        for (int i = 0;i < len; i++) bit[i] &= b.bit[i];
    }
    
    void operator |= (const Bit &b)  {
        for (int i = 0;i < len; i++) bit[i] |= b.bit[i];
    }
    
    void operator ^= (const Bit &b)  {
        for (int i = 0;i < len; i++) bit[i] ^= b.bit[i];
    }
    
    Bit operator << (const int t) const {
        Bit res; int high = t >> 6, low = t & 63;
        ull last = 0;
        for (int i = 0;i + high < len; i++) {
            res.bit[i + high] = (last | (bit[i] << low));
            if (low) last = (bit[i] >> (64 - low));
        }
        return res;
    }
    
    Bit operator >> (const int t) const {
        Bit res; int high = t >> 6, low = t & 63;
        ull last = 0;
        for (int i = len - 1;i >= high; i--) {
            res.bit[i - high] = last | (bit[i] >> low);
            if (low) last = bit[i] << (64 - low);
        }
        return res;
    }
    
    void operator <<= (const int t) {
        int high = t >> 6, low = t & 63;
        for (int i = len - high - 1; ~i; i--) {
            bit[i + high] = (bit[i] << low);
            if (low && i) bit[i + high] |= bit[i - 1] >> (64 - low);
        }
        for (int i = 0;i < high; i++) bit[i] = 0;
    }
};