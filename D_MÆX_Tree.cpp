
/*
" If you set your goals ridiculously high, even if its a failure
You will fail above everyone else's success" 
*/

#pragma GCC optimize("O1")
#pragma GCC optimize("O2")
#pragma GCC optimize("Ofast")
#pragma GCC target("avx,avx2,fma")
#pragma GCC optimize("O3,unroll-loops")
#include <bits/stdc++.h>
using namespace std;

#include<ext/pb_ds/assoc_container.hpp>
#include<ext/pb_ds/tree_policy.hpp>
using namespace __gnu_pbds;
using namespace __gnu_cxx;
typedef tree<int,null_type,less<int>,rb_tree_tag,tree_order_statistics_node_update> ind_set ;
typedef tree<int,null_type,less_equal<int>,rb_tree_tag,tree_order_statistics_node_update> ind_mset ;
typedef map < int , tree < int , null_type , less < int > , rb_tree_tag , tree_order_statistics_node_update  > > ind_map;

typedef unsigned long long ull;
typedef long double ld;
typedef long long ll;
typedef pair<int, int> pii;
typedef vector<int> vi;
typedef vector<vi> vii;
typedef vector<pii> vpii;
typedef pair<ll, ll> pll;
typedef vector<ll> vl;
typedef vector<pll> vpll;
typedef vector<vl> vll;
template<class T> using maxpq = priority_queue<T>;
template<class T> using minpq = priority_queue<T, vector<T>, greater<T>>;
template<class T> using vec1 = vector<T>;
template<class T> using vec2 = vector<vector<T>>;
template<class T> using vec3 = vector<vector<vector<T>>>;

#define FOR(i, a, b) for (int i = (a); i < (b); i++)
#define REP(i,b,a)  for(int i =((b)-1);i>=(a);i--)
#define mem(a,x)   memset(a,x,sizeof(a))
#define pb push_back
#define pf push_front
#define ff first
#define ss second
#define sz(x) ((int)(x).size())
#define all(v) v.begin(), v.end()
#define rall(v) (v).rbegin(),(v).rend()
#define eb emplace_back
#define mp make_pair
#define mt make_tuple
#define fbo(a) find_by_order(a) //will give a-th largest element
#define ook(a) order_of_key(a) //will give no. of elements strictly lesser than a
#define sz(x) ((int)(x).size())
#define nzl(x) __builtin_clzll(x)
#define nzr(x) __builtin_ctzll(x)
#define setbits(x) __builtin_popcountll(x)
#define setbitsParity(x) __builtin_parityll(x) // 1 -> odd else 0 if even
#define umap unordered_map
#define uset unordered_set
#define nl "\n"
#define PI atan(1)*4
#define E 2.71828
#define yes {cout << "Yes" << endl; return;}
#define no {cout << "No" << endl; return;}
#define YES {cout << "Yes" << endl;}
#define NO {cout << "No" << endl;}
#define nyet {cout<<"-1"<<endl;return;}
#define mxe(v)  (*max_element(v.begin(),v.end()))
#define mne(v)  (*min_element(v.begin(),v.end()))
#define unq(v)  v.resize(distance(v.begin(), unique(v.begin(), v.end())));
#define ub upper_bound
#define lb lower_bound
#define LB(c, x) distance((c).begin(), lower_bound(all(c), (x)))
#define UB(c, x) distance((c).begin(), upper_bound(all(c), (x)))
#define UNIQUE(x) \
  sort(all(x)), x.erase(unique(all(x)), x.end()), x.shrink_to_fit()
#define outt(a) \
      FOR(i,1,sz(a))            \
      cout << a[i] << " "; \
      cout << endl;
#define inn(a) \
      FOR(i,1,sz(a))            \
      cin>>a[i];
#define FAST_AF_BOI                \
    ios_base ::sync_with_stdio(0); \
    cin.tie(0);               \
    cout.tie(0);

struct custom_hash {
  static uint64_t splitmix64(uint64_t x) {
      x += 0x9e3779b97f4a7c15;//abk
      x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9;
      x = (x ^ (x >> 27)) * 0x94d049bb133111eb;
      return x ^ (x >> 31);
  }
  size_t operator()(uint64_t x) const {
      static const uint64_t FIXED_RANDOM = chrono::steady_clock::now().time_since_epoch().count();
      return splitmix64(x + FIXED_RANDOM);
  }
};
typedef gp_hash_table<long long,long long,custom_hash> fast_map;
typedef unordered_map<long long,long long,custom_hash> safe_map;

// ================================== i/o module ==================================
template <class T> void _print(T x){cerr<<x;};
template <class T, class V> void _print(pair <T, V> p);
template <class T> void _print(vector <T> v);
template <class T> void _print(set <T> v);
template <class T, class V> void _print(map <T, V> v);
template <class T> void _print(multiset <T> v);
template <class T, class V> void _print(pair <T, V> p) {cerr << "{"; _print(p.ff); cerr << ","; _print(p.ss); cerr << "}";}
template <class T> void _print(vector <T> v) {cerr << "[ "; for (T i : v) {_print(i); cerr << " ";} cerr << "]";}
template <class T> void _print(set <T> v) {cerr << "[ "; for (T i : v) {_print(i); cerr << " ";} cerr << "]";}
template <class T> void _print(multiset <T> v) {cerr << "[ "; for (T i : v) {_print(i); cerr << " ";} cerr << "]";}
template <class T, class V> void _print(map <T, V> v) {cerr << "[ "; for (auto i : v) {_print(i); cerr << " ";} cerr << "]";}
template<class T>void read(T &x){x=0;int f=0;char ch=getchar();while(ch<'0' || ch>'9')f|=(ch=='-'),ch=getchar();while(ch>='0' && ch<='9')x=(x<<3)+(x<<1)+(ch^48),ch=getchar();x=f? -x:x;return ;}
template<typename typC,typename typD> istream &operator>>(istream &cin,pair<typC,typD> &a) { return cin>>a.first>>a.second; }
template<typename typC> istream &operator>>(istream &cin,vector<typC> &a) { for (auto &x:a) cin>>x; return cin; }
template<typename typC,typename typD> ostream &operator<<(ostream &cout,const pair<typC,typD> &a) { return cout<<a.first<<' '<<a.second; }
template<typename typC,typename typD> ostream &operator<<(ostream &cout,const vector<pair<typC,typD>> &a) { for (auto &x:a) cout<<x<<'\n'; return cout; }
template<typename typC> ostream &operator<<(ostream &cout,const vector<typC> &a) { int n=a.size(); if (!n) return cout; cout<<a[0]; for (int i=1; i<n; i++) cout<<' '<<a[i]; return cout; }
template <class T> inline vector<T>& operator--(vector<T>& v) { for(auto &x:v) --x; return v; }
template <class T> inline vector<T>& operator++(vector<T>& v) { for(auto &x:v) ++x; return v; }
template <class T> inline vector<T>& operator^=(vector<T>& v,int y) { for(auto &x:v) x^=y; return v; }
inline string& operator^=(string& s,int y) { for(auto &x:s)x=((x-'0')^y)+'0' ; return s; }

void dgb_out () { cerr << endl; }
template < typename Head, typename... Tail >
void dgb_out ( Head H, Tail... T) { cerr <<' ' << H; dgb_out (T...); }
#ifndef ONLINE_JUDGE
#define dbg(...) cerr << "(" << #__VA_ARGS__ << "):", dgb_out(__VA_ARGS__) 
#else
#define dbg(...)
#endif

// ================================================================================

//`````````````````````````````````````````````IMP FUNCTIONS``````````````````````````````````````````````````````
ll ceil(ll a,ll b){return (a+b-1)/b;}
int log_2(ull i){return i?nzl(1)-nzl(i):-1;}
ll gcd(ll a, ll b) {if (b > a) {return gcd(b, a);} if (b == 0) {return a;} return gcd(b, a % b);}
ll bin_expo(ll a, ll b, ll mod) {ll res = 1;a%=mod;if(a==0)return 0;while (b > 0) {if (b & 1)res = (res * a) % mod; a = (a * a) % mod; b = b >> 1;} return res;}
ll bin_mul(ll a, ll b, ll mod) {ll res = 0; while (b > 0) {if (b & 1)res = (res + a) % mod; a = (a + a) % mod; b = b >> 1;} return res;}
void extendgcd(ll a, ll b, ll*v) {if (b == 0) {v[0] = 1; v[1] = 0; v[2] = a; return ;} extendgcd(b, a % b, v); ll x = v[1]; v[1] = v[0] - v[1] * (a / b); v[0] = x; return;} //pass an arry of size1 3
ll mminv(ll a, ll b) {ll arr[3]; extendgcd(a, b, arr); return arr[0];} //for non prime b
ll mminvprime(ll a, ll b) {return bin_expo(a, b - 2, b);}
ll ncr(ll n, ll r, ll m, ll *fact, ll *ifact) {if(n<r)return 0;ll val1 = fact[n]; ll val2 = ifact[n - r]; ll val3 = ifact[r];if(n < r) return 0;if(n == r || r == 0) return 1;if(r<0) return 0; return (((val1 * val2) % m) * val3) % m;}
void google(int t) {cout << "Case #" << t << ": ";}
vl sieve(int n) {int*arr = new int[n + 1](); vl vect; for (int i = 2; i <= n; i++)if (arr[i] == 0) {vect.push_back(i); for (int j = 2 * i; j <= n; j += i)arr[j] = 1;} return vect;}
ll mod_add(ll a, ll b, ll m) {a = a % m; b = b % m; return (((a + b) % m) + m) % m;}
ll mod_mul(ll a, ll b, ll m) {a = a % m; b = b % m; return (((a * b) % m) + m) % m;}
ll mod_sub(ll a, ll b, ll m) {a = a % m; b = b % m; return (((a - b) % m) + m) % m;}
ll mod_div(ll a, ll b, ll m) {a = a % m; b = b % m; return (mod_mul(a, mminvprime(b, m), m) + m) % m;}  //only for prime m
ll phin(ll n) {ll number = n; if (n % 2 == 0) {number /= 2; while (n % 2 == 0) n /= 2;} for (ll i = 3; i*i <= n; i += 2) {if (n % i == 0) {while (n % i == 0)n /= i; number = (number / i * (i - 1));}} if (n > 1)number = (number / n * (n - 1)) ; return number;} //O(sqrt(N))
ll large_expo(ll a,ll b,ll c,ll m) {return bin_expo(a,bin_expo(b,c,phin(m)),m);} //(a^b^c)%M 
ll large_expo_prime(ll a,ll b,ll c,ll m) {return bin_expo(a,bin_expo(b,c,m-1),m);} //(a^b^c)%M when m is prime
template<class T>vector<T> prefixSum(vector<T> v, bool flag){vector<T> ans;T sum = 0;if (flag){for (auto &e : v){sum += e;ans.eb(sum);}}else{ans.pb(0);REP(i, v.size(), 0){sum += v[i];ans.eb(sum);}reverse(all(ans));}return ans;}
ll ffs(ll n){if(n==0)return -1;return log2(n & -n);}
template<class T> bool ckmin(T& a, const T& b) { return b < a ? a = b, 1 : 0; }
template<class T> bool ckmax(T& a, const T& b) { return a < b ? a = b, 1 : 0; }
mt19937_64 rng(chrono::steady_clock::now().time_since_epoch().count());
ll uid(ll l, ll r) {return uniform_int_distribution<ll>(l, r)(rng);} 
//````````````````````````````````````````````````````````````````````````````````````````````````````````````

const ll INF=1e18+10;
const ll inf=1e9+10;
const ll N=1e5+10;
const ll mod1=(1e9+7);
const ll mod2=(998244353);

/*
    n    - 1e5
    q    - 1e5
    a[i] - 2e5
    TC - O( nlogn + q * ( n^(2/3) + sqrt(n) )) - O( 2.5e8 )
*/

// #define int long long

const int L=17;

int SZ,SZ_t,SZ_n;
int n,k,timer,m,t,q;
vec2<int>adj(N);
int a[N],tmp[N],ans[N],cnt[N];
int tin[N],tout[N],dpth[N];
int freq[N],ifreq[N],ok[N];
int lift[N][L];
vector<int>val;
uset<int,custom_hash>st;

void dfs(int u,int p){
    val.pb(u);
    tin[u]=++timer;
    lift[u][0]=p;
    for(auto &v:adj[u]){
        if(v==p)continue;
        dpth[v]=dpth[u]+1;
        dfs(v,u);
    }
    val.pb(u);
    tout[u]=++timer;
}

void Compute(){
    val.pb(0);
    dfs(1,0);
    FOR(j,1,L)
        FOR(i,1,n+1)
            lift[i][j]=lift[lift[i][j-1]][j-1];
}

int kthanc(int x,int d){
    FOR(i,0,L)
        if(d>>i&1)x=lift[x][i];
    return x;
}

int lca(int u,int v){
    if(dpth[u]>dpth[v])swap(u,v);
    v=kthanc(v,dpth[v]-dpth[u]);
    if(u==v)return u;
    REP(i,L,0)
        if(lift[v][i]!=lift[u][i])u=lift[u][i],v=lift[v][i];
    return lift[u][0];
}
int64_t peanoOrder(unsigned x, unsigned y, unsigned z, int m) {
    static int a[3][13];
    memset(a, 0, sizeof(a));
    int sum0 = 0, sum1 = 0, sum2 = 0;
    int ptr = m - 1;
    while (x) {
        a[0][ptr] = x % 3;
        sum0 += a[0][ptr];
        ptr--;
        x /= 3;
    }
    ptr = m - 1;
    while (y) {
        a[1][ptr] = y % 3;
        sum1 += a[1][ptr];
        ptr--;
        y /= 3;
    }
    ptr = m - 1;
    while (z) {
        a[2][ptr] = z % 3;
        sum2 += a[2][ptr];
        ptr--;
        z /= 3;
    }
 
    for (int i = m - 1; i >= 0; i--) {
        sum2 -= a[2][i];
        if ((sum0 + sum1) & 1) a[2][i] = 2 - a[2][i];
 
        sum1 -= a[1][i];
        if ((sum0 + sum2) & 1) a[1][i] = 2 - a[1][i];
 
        sum0 -= a[0][i];
        if ((sum1 + sum2) & 1) a[0][i] = 2 - a[0][i];
    }
 
    int64_t num = 0, base = 1;
    for (int j = m - 1; j >= 0; j--) {
        num += base * a[2][j];
        base *= 3;
        num += base * a[1][j];
        base *= 3;
        num += base * a[0][j];
        base *= 3;
    }
    return num;
}
struct que{
    int t,l,r,i,node,k;
    ll ord;
    que(){}
    que(int t,int l,int r,int i,int k,int node=0){
        this->t=t;
        this->l=l;
        this->r=r;
        this->i=i;
        this->k=k;
        this->node=node;
        this->ord=peanoOrder(l,r,t,13);
    }
    bool operator < (const que &u) const {
        return ord<u.ord;
    }
}Q[N];

struct upd{
    int node,prev,cur;
    upd(){}
    upd(int node,int prev,int cur){
        this->node=node;
        this->prev=prev;
        this->cur=cur;
    }
}U[N];

void add(int x){
    ifreq[freq[x]]--;
    if(!ifreq[freq[x]])st.erase(freq[x]);
    freq[x]++;
    if(!ifreq[freq[x]])st.insert(freq[x]);
    ifreq[freq[x]]++;
}

void del(int x){
    ifreq[freq[x]]--;
    if(!ifreq[freq[x]])st.erase(freq[x]);
    freq[x]--;
    if(!ifreq[freq[x]] and freq[x])st.insert(freq[x]);
    ifreq[freq[x]]++;
}

void insert(int node){
    if(ok[node]){
        del(a[node]);
        ok[node]=0;
    }
    else {
        add(a[node]);
        ok[node]=1;
    }
}

void erase(int node){
    if(ok[node]){
        del(a[node]);
        ok[node]=0;
    }
    else {
        add(a[node]);
        ok[node]=1;
    }
}

void ahead(int t,int l,int r){
    int node=U[t].node,prev=U[t].prev,cur=U[t].cur;
    int in=tin[node],out=tout[node];
    a[node]=cur;
    int ok1=in>=l and in<=r,ok2=out>=l and out<=r;
    if(ok1^ok2){
        del(prev);
        add(cur);
    }
}

void behind(int t,int l,int r){
    int node=U[t].node,prev=U[t].prev,cur=U[t].cur;
    int in=tin[node],out=tout[node];
    a[node]=prev;
    int ok1=in>=l and in<=r,ok2=out>=l and out<=r;
    if(ok1^ok2){
        del(cur);
        add(prev);
    }
}
void transcendent(int tc)
{
    cin>>n;
    FOR(i,1,n+1)cin>>a[i],tmp[i]=a[i];
    FOR(_,0,n-1){
        int u,v;cin>>u>>v;
        adj[u].pb(v);
        adj[v].pb(u);
    }
    Compute();
    cin>>m;
    FOR(_,1,m+1){
        int type;cin>>type;
        if(type==1){
            int u,x;cin>>u>>x;
            U[++t]=upd(u,tmp[u],x);
            tmp[u]=x;
        }else{
            int u,v;cin>>u>>v>>k;
            if(tin[u]>tin[v])swap(u,v);
            int l=lca(u,v);
            if(l==u){
                // [tin[u],tin[v]]
                q++;
                Q[q]=que(t,tin[u],tin[v],q,k);
            }else{
                // [tout[u],tin[v]] + lca
                q++;
                Q[q]=que(t,tout[u],tin[v],q,k,l);
            }
        }
    }
    SZ=sqrtl(n);
    while(SZ*(SZ+1)/2<=n)SZ++;
    sort(Q+1,Q+q+1);
    for(int l=1,r=0,t=0,i=1;i<=q;i++){
        while(r<Q[i].r)insert(val[++r]);
        while(l>Q[i].l)insert(val[--l]);
        while(r>Q[i].r)erase(val[r--]);
        while(l<Q[i].l)erase(val[l++]);
        while(t<Q[i].t)ahead(++t,l,r);
        while(t>Q[i].t)behind(t--,l,r);
        if(Q[i].node)add(a[Q[i].node]);
        assert(sz(st)<=SZ);
        k=Q[i].k;
        cnt[0]=1;
        for(auto &x:st){
            assert(ifreq[x]);
            cnt[x%k]+=ifreq[x];
        }
        ll mx=INF;
        FOR(i,0,min(SZ+1,k)){
            ckmin(mx,i+cnt[i]*1ll*k);
        }
        ans[Q[i].i]=mx;
        cnt[0]=0;
        for(auto &x:st){
            cnt[x%k]=0;
        }
        if(Q[i].node)del(a[Q[i].node]);
    }
    FOR(i,1,q+1)cout<<ans[i]<<nl;
}

static void read(){
    freopen("input.txt","r",stdin);
    freopen("output.txt","w",stdout);
}

int32_t main()
{
    //read();
    FAST_AF_BOI
    auto begin = std::chrono::high_resolution_clock::now();
    cout << fixed << setprecision(0);
    cerr << fixed << setprecision(0);
    //clock_t timer;
    //timer = clock();
    //PreComp();
    int test=1;
    // cin >> test;
    FOR(tc,1,test+1)
    {
        // cerr<<endl<<"----Test:"<<tc<<"----"<<endl;
        transcendent(tc);
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    //cerr << "Time measured: " << elapsed.count() * 1e-9 << " seconds.\n"; 
    return 0;
}
