#include<bits/stdc++.h>
using namespace std;
#define int long long
 
signed main(){
    ios_base::sync_with_stdio(false);cin.tie(NULL);cout.tie(NULL);
    int n,k;
    cin >> n >> k;
    vector<int> a(n);
    for(int i = 0; i < n; i++){
        cin >> a[i];
    }
    int idx = -1;
    for(int i = 0; i < n; i++){
        if(a[i] == k){
            idx = i;
            break;
        }
    }
    int cnt = 1;
    for(int i = idx + 1; i < n; i++){
        if(a[i] > k) cnt++;
        else break;
    }
    for(int i = idx - 1; i >= 0; i--){
        if(a[i] > k) cnt++;
        else break;
    }
    cout << n - cnt << endl;
  return 0;
}
