// by Siddhid Saha (2112010)


#include <bits/stdc++.h>
using namespace std;

#define INF 1e18
#define endl "\n"
#define pb push_back
#define ppb pop_back
#define mp make_pair
#define PI atan(1)*4
#define set_bits __builtin_popcountllO
#define all(x) (x).begin(), (x).end()
#define vi vector<int>
#define vll vector<ll>
#define pll pair<ll,ll>
#define rvsort(a) sort(all(a),greater<int>())
#define read(a,n) for(int i = 0 ; i < n ; i ++){ cin >> a[i];}
#define printv(a) for(auto it: a){cout << it << " ";} cout << endl;
#define ms(arr, v) memset(arr, v, sizeof(arr))

typedef long long ll;
typedef unsigned long long ull;
typedef long double lld;

mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());

#ifndef ONLINE_JUDGE
#include "/Users/templates/debug.h"
#else
#define dbg(x...)
#endif
/*---------------------------------------------------------------------------------------------------------------------------*/
ll gcd(ll a, ll b) {if (b > a) {return gcd(b, a);} if (b == 0) {return a;} return gcd(b, a % b);}
ll expo(ll a, ll b, ll mod) {ll res = 1; while (b > 0) {if (b & 1)res = (res * a) % mod; a = (a * a) % mod; b = b >> 1;} return res;}
void extendgcd(ll a, ll b, ll*v) {if (b == 0) {v[0] = 1; v[1] = 0; v[2] = a; return ;} extendgcd(b, a % b, v); ll x = v[1]; v[1] = v[0] - v[1] * (a / b); v[0] = x; return;} //pass an arry of size1 3
ll mminv(ll a, ll b) {ll arr[3]; extendgcd(a, b, arr); return arr[0];} //for non prime b
ll mminvprime(ll a, ll b) {return expo(a, b - 2, b);}
void google(int t) {cout << "Case #" << t << ": ";}
vector<ll> sieve(int n) {int*arr = new int[n + 1](); vector<ll> vect; for (int i = 2; i <= n; i++)if (arr[i] == 0) {vect.push_back(i); for (int j = 2 * i; j <= n; j += i)arr[j] = 1;} return vect;}
ll phin(ll n) {ll number = n; if (n % 2 == 0) {number /= 2; while (n % 2 == 0) n /= 2;} for (ll i = 3; i <= sqrt(n); i += 2) {if (n % i == 0) {while (n % i == 0)n /= i; number = (number / i * (i - 1));}} if (n > 1)number = (number / n * (n - 1)) ; return number;} //O(sqrt(N))
ll uid(ll l, ll r) {return uniform_int_distribution<ll>(l, r)(rng);} 
/*--------------------------------------------------------------------------------------------------------------------------*/
//const int mod = 1e9 + 7;
//const int mod = 998244353; 



int n, l;
vector<vector<int>> adj;
vll value;
struct LCA{
    vector<vector<int>> g;
    vector<vector<int>> up;
    vll tin, tout, depth, sumv;
    int timer = 0;
    inline void init(int n){
        n++;
        g.resize(n);
        tin.resize(n);
        tout.resize(n);
        depth.resize(n);
        sumv.resize(n);
        up.assign(n, vector<int>(31, -1));
    }
    LCA(){};
    LCA(int n){
        init(n);
    }

    inline void add_edge(int u, int v){
        g[u].push_back(v);
        g[v].push_back(u);
    }

    inline void dfs(int v, int p , int dis, ll val){
        tin[v] = ++timer;
        depth[v] = dis;
        val += value[v];
        sumv[v] = val; 
        up[v][0] = p;
        for (int i = 1; i < 31; i++){
            up[v][i] = up[up[v][i - 1]][i - 1];
        }

        for (int u : g[v]){
            if (u != p){
                dfs(u, v , dis+1, val);
            }
        }
        tout[v] = ++timer;
    }

    void dfs(int root = 1){
        dfs(root, root , 0, 0);
    }

    inline bool upper(int u, int v){
        return tin[u] <= tin[v] && tout[v] <= tout[u];
    }

    int lca(int u, int v){
        if (upper(u, v)) return u;
        if (upper(v, u)) return v;
        for (int i = 30 ; i >= 0; i--){
            if (!upper(up[u][i], v)) u = up[u][i];
        }
        return up[u][0];
    }

    inline int dist(int u, int v){
        return depth[u] + depth[v] - 2 * depth[lca(u, v)];
    }
    ll valuv(ll u, ll v){
        return sumv[u] + sumv[v] - 2*sumv[lca(u, v)]+value[lca(u, v)];
    }

};

void solve()
{
    cin >> n;
    assert(n >= 1 && n <= 1e5);
    adj.resize(n+1);
    LCA lc(n);
    for(int i = 0 ; i < n-1 ; i ++){
        ll u, v; cin >> u >> v;
        assert(u != v && u >= 1 && u <= n && v >= 1 && v <= n);
        lc.add_edge(u, v);
    }
    value.resize(n+1);
    for(int i = 1 ; i<= n ; i ++) cin >> value[i];
    for(int i = 1 ; i <= n ; i ++){
        assert(value[i] >= 0 && value[i] <= 1e9);
    }
    lc.dfs();
    ll q; cin >> q;
    assert(q >= 1 && q <= 1e5);
    while(q--){
        ll u, v; cin >> u >> v;
        assert(u != v && u >= 1 && u <= n && v >= 1 && v <= n);
        ll lca = lc.lca(u ,v);
        ll d = lc.dist(u, v), du = lc.dist(u, lca);
        ll l = 0, r = d;
        while(l < r){
            ll mid = (l+r)/2;
            auto check = [&](ll mid) mutable {
                if(mid > du){
                    ll x = v;
                    ll up = d-mid-1;
                    for(int i = 30 ; i >= 0 ; i--){
                        if(up&(1LL<<i)){
                            x = lc.up[x][i];
                        }
                    }
                    ll vv = lc.valuv(v, x);
                    ll totv = lc.valuv(u, v);
                    ll remv = totv-vv;
                    ll totnum = d+1, numv = mid+1;
                    ll remnum = totnum-numv;
                    vv += ((remnum)*(remnum+1))/2;
                    remv += ((numv)*(numv+1))/2;
                    return remv > vv;
                }else{
                    ll x = u;
                    ll up = mid;
                    for(int i = 30 ; i >= 0 ; i--){
                        if(up&(1LL<<i)){
                            x = lc.up[x][i];
                        }
                    }
                    ll vv = lc.valuv(u, x);
                    ll totv = lc.valuv(u, v);
                    ll remv = totv-vv;
                    // dbg(u, x);
                    // dbg(vv, remv);
                    ll totnum = d+1, numv = mid+1;
                    ll remnum = totnum-numv;
                    vv += ((numv)*(numv+1))/2;
                    remv += ((remnum)*(remnum+1))/2;
                    return vv > remv;
                }
            };
            if(check(mid)){
                r = mid;
            }else{
                l = mid+1;
            }
        }
        // dbg(d, l);
        ll hh = d-l+1;
        cout << hh << endl;
    }



}

int main() {

ios::sync_with_stdio(0);
cin.tie(0);
cout.tie(0);


ll t = 1;
// cin >> t;
for(int i = 1 ; i <= t ; i++){
//google(i);
solve();
}
return 0;
}








