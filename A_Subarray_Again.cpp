//Author : RAJ ROY
//scholar id : 2212002

#include <bits/stdc++.h>
using namespace std;

// #define endl "\n"
#define ll long long
vector<long long >result;


#pragma GCC optimize("O1")
#pragma GCC optimize("O2")
#pragma GCC optimize("Ofast")
#pragma GCC target("avx,avx2,fma")
#pragma GCC optimize("O3,unroll-loops")
#pragma comment(linker,"/stack:200000000")
#pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx,avx,tune=native")
#pragma GCC optimize("Ofast,no-stack-protector,no-stack-protector,fast-math")
/*****************************************************************************************************************************************/


void solve() {
    int n, k;
    cin >> n >> k;
    vector<int> arr(n+1);
    for (int i = 1; i <= n; ++i) {
        cin >> arr[i];
    }
    vector<vector<int>> pre(n+1, vector<int>(32, 0));
    for (int i = 1; i <= n; ++i) {
        for (int j = 0; j < 32; ++j) {
            if ((1 << j) & arr[i]) {
                pre[i][j] = 1;
            }
            pre[i][j] += pre[i - 1][j];
        }
    }

    long long ans = 0;
    for (int i = 1; i <= n; ++i) {
        int l = i, r = n;
        int an = -1;
        while (l <= r) {
            int mid = (l + r) / 2;
            int d = 0;
            for (int j = 0; j < 32; ++j) {
                if (pre[mid][j] - pre[i-1][j] > 0) {
                    d |= (1 << j);
                }
            }
            if (d >= k) {
                if (d == k) an = mid;
                r = mid - 1;
            } else {
                l = mid + 1;
            }
        }

        if (an != -1) {
            l = an, r = n;
            int an1 = an;
            while (l <= r) {
                int mid = (l + r) / 2;
                int d = 0;
                for (int j = 0; j < 32; ++j) {
                    if (pre[mid][j] - pre[i-1][j] > 0) {
                        d |= (1 << j);
                    }
                }
                if (d > k) {
                    r = mid - 1;
                } else {
                    if (d == k) an1 = mid;
                    l = mid + 1;
                }
            }
            ans += (an1 - an + 1);
        }
    }
    cout << ans << endl;

}
int main(){
 
    #ifndef ONLINE_JUDGE
    freopen("input.txt", "r", stdin);
    freopen("output.txt", "w", stdout);
    #endif
 
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int t = 1;
    //int t; cin >> t;
    for(int i = 1; i <= t; i++){
        //cout << "Case #" << i << ":"<< " ";
        solve();
    }

    return 0;
}
