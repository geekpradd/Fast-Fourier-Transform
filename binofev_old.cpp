#include <bits/stdc++.h>
#define int long long
#define MOD 998244353
using namespace std;

int power(int a, int b){
	if (b==0)
		return 1;
	int num = power(a, b/2);
	num = (num*num)%MOD;
	if (b%2)
		return (num*a)%MOD;
	else
		return num;
}

int inverse(int a){
	return power(a, MOD-2);
}

void get_coefficients(vector<int> &poly, int l, int r){
	if (l==r){
		poly.resize(2);
		poly[0] = -1*(l); poly[1] = 1;
		if (poly[0]!=0)
			poly[0] = poly[0] + MOD;
		return;
	}
	int mid = (l+r)/2;
	vector<int> left(mid - l + 2, 0);
	vector<int> right(r - (mid + 1) + 2, 0);
	get_coefficients(left, l, mid);
	get_coefficients(right, mid+1, r);
	// Below we need to multiply left with right to get the final polynomial
	// naive implementation will be O(r^2) 
	// using FFT with NTT we will bring it down to O(rlogr)
	// using naive now

	for (int i=0; i<left.size(); ++i){
		for (int j=0; j<right.size(); ++j){
			int term = (left[i]*right[j])%MOD;
			poly[i+j]= (poly[i+j] + term)%MOD;
		}
	}


}

signed main(){
	cin.tie(NULL); ios_base::sync_with_stdio(false);
	int t; cin >> t;
	// CASE P = 1 needs to be handled differently 
	while (t--){
		int n, p, r; cin >> n >> p >> r;
		// if (p==1){
		// 	if (r==1)
		// 		cout << (n)%MOD << endl;
		// 	else
		// 		cout << 0 << endl;

		// 	continue;
		// }
		int fact = 1;
		for (int i=1; i<=r; ++i)
			fact = (fact*i)%MOD;
		fact = inverse(fact);

		vector<int> poly(r+1, 0);
		get_coefficients(poly, 0, r-1);

		int ans = 0;
		for (int i=1; i<=r; ++i){
			// cout << "poly " << i << " is " << poly[i] << endl;
			int num = (power(p, (n+1)*i) - 1)%MOD;
			int denom = (power(p, i) - 1)%MOD;
			if (denom != 0){
				denom = inverse(denom);
				ans = ( ans + ((num*denom)%MOD*poly[i])%MOD)%MOD;
			}
			else {
				// common ratio of GP is congruent to 1
				ans = (ans + ((n+1)%MOD*poly[i])%MOD)%MOD;
			}
		}
		ans = (ans*fact)%MOD;
		cout << ans << endl;
	}

	return 0;
}