#include<iostream>
 
using namespace std;
#define endl '\n'
#define pi acos(-1)
#define euler_e 2.71828
#define pii pair<ll int ,ll int>
#define pb push_back
#define readv(v) for(auto &x: v) cin>>x  
#define ll long long
//#define all(str) str.begin(), str.end()
#define declear(x) decltype(x)
#define ull unsigned long long
#define io freopen("input.txt","r",stdin);freopen("output.txt","w",stdout);
#define read freopen("input.txt","r",stdin);
#define boost ios_base::sync_with_stdio(false);cin.tie(0);cout.tie(0);
#define errd 0.000000001
//int kx[] = {+2, +1, -1, -2, -2, -1, +1, +2};
//int ky[] = {+1, +2, +2, +1, -1, -2, -2, -1};
int fx[] = {+0, +0, +1, -1, -1, +1, -1, +1};
int fy[] = {-1, +1, +0, +0, +1, +1, -1, -1};
const ll int infLL=9223372036854775103;
const int inf=2147483643;
ll int mod=(ll)(1e9+7);
ll int ceil(ll int a,ll int b){return (a+b-1)/b;}
ll int min(ll int a, ll int b){if(a>b)return b;else return a;}
bool bit_check(ll int a,int i){
  if((a & (1LL<<i)))return 1;
  return 0;
}
ll int bit_toggle(ll int a,int i){
  return (a^(1LL<<i));
}
ll int bit_sum(ll int a,int i){
  return a+(1LL<<i);
}
ll int bit_sub(ll int a,int i){
  return a-(1LL<<i);
}
ll int mod_power(ll int x,ll int y){//x^y%p
  ll int p=mod;
  ll int res = 1;
  x = x % p;
  while (y > 0) { 
    if (y & 1)res = (res*x) % p;
    y = y>>1;x = (x*x) % p;
  }
  return res;
}
ll int power_of(ll int a,int b){
  if(a==0)return -1;
  return 1+power_of(a/b,b);
}
ll power(ll int a, ll int b) {
    if(a==1)return 1;
    ll int res = 1;
    while (b > 0) {
        if (b & 1)
            res = res * a;
        a = a * a;
        b >>= 1;
    }
    return res;
}
 
//auto cmp = [](pii left, pii right) { return cmp2(left,right); };
void mycode();
int main(){
    boost;
    //io;
    //freopen("output.txt","w",stdout);
    //read;
    mycode();
    //cerr<<"Ended\n";
    return 0;
}
#include <vector>
#include <set>
#include <algorithm>
#include <cstring>
template <typename T>
int intersec_size(const set<T>& a, const set<T>& b){
    set<T> result = a;
    result.insert(b.begin(), b.end());
    return a.size()+b.size()-result.size();
}
template <typename T>
int union_size(const set<T>& a, const set<T>& b){
    set<T> result = a;
    result.insert(b.begin(), b.end());
    return result.size();
}
 
int arr2[40][3];
int visited[40][3];
vector<set<pair<int,int>>> vc;
set<pair<int,int>> tmp;
 
int n, m, k;
int dfs(int x, int y){
    if(arr2[x][y] == 1)return -90000;
    if(visited[x][y] == 1)return 0;
    if(arr2[x][y] == -1){
        return 0;
    }
    else{
        visited[x][y] = 1;
        int sum = 0;
        for(int j=0;j<4;j++){
            if((x+fx[j])>=0 && (x+fx[j])<n && y+fy[j]>=0 && y+fy[j]<2){
                sum+=dfs(x+fx[j],y+fy[j]);
            }
        }
        return 1+sum;
    }
}
 
void mycode(){
    cin>>n>>m>>k;
    int a, b;
    int arr[n+2][3]={0};
    memset(arr, 0, sizeof(arr));
    for(int i=0;i<m;i++){
        cin>>a>>b;
        arr[b-1][a-1] = 1;
    }
    set<pair<int,int>>st;
    for(int i=0;i<n;i++){
        for(int l = 0;l<2;l++){
            for(int j=0;j<4;j++){
                if(arr[i][l] && (i+fx[j])>=0 && (i+fx[j])<n && l+fy[j]>=0 && l+fy[j]<2 && !arr[i+fx[j]][l+fy[j]]){
                    st.insert({i+fx[j],l+fy[j]});
                }
            }
        }
    }
    for(int i=0;i<n;i++){
        arr2[i][0]=arr[i][0];
        arr2[i][1]=arr[i][1];
        visited[i][0] = 0;
        visited[i][1] = 0;
    }
    int at = 0;
    pii positions[100];
    for(int i=0;i<n;i++){
        if(arr[i][0]==0)positions[at++] = {i,0};
        if(arr[i][1]==0)positions[at++] = {i,1};
    }

    vector<int> permo(at);
    for(int i=0;i<k;i++)permo[i]=-1;
    for(int i=k;i<at;i++)permo[i]=0;
    sort(permo.begin(), permo.end());
    int ansmx = 0;
    do {
        memset(visited, 0, sizeof(visited));
        at = 0;
        int ans2=0, ans0, ans1;
        for(auto &x:permo){
            arr2[positions[at].first][positions[at++].second] = x;
        }
        for(int i=0;i<n;i++){
            if(arr2[i][0] == 0 && !visited[i][0]){
                ans0 = dfs(i,0);
                if(ans0>0)ans2+=ans0;
            }
            if(arr2[i][1] == 0 && !visited[i][1]){
                ans1 = dfs(i,1);
                if(ans1>0)ans2+=ans1;
            }
        }
        ansmx = max(ansmx, ans2);
    } while (next_permutation(permo.begin(), permo.end()));
    cout<<min(ansmx+k,2*n-m)<<endl;
}
