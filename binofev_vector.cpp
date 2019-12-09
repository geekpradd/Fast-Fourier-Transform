#include <bits/stdc++.h>
#define int long long
#define primitive 3
#define primitive_i 332748118
#define MOD 998244353
using namespace std;

int power(int a, int b){
	if (b==0)
		return 1;
	int num = power(a, b/2);
	num = (int)((1LL*num*num)%MOD);
	if (b%2)
		return (int)((1LL*num*a)%MOD);
	else
		return num;
}
int inverse(int a){
	return power(a, MOD-2);
}

// setting up constants
// int primitive = 3; // primitive root of MOD, existence proved by Euler
// int root = power(primitive, 7*17); // root will be the 2^23 root of unity
// // in the ring mod MOD
// int root_i = inverse(root);
// int limit = 1 << 23;


// goal: NTT implementation of fast fourier transform 
void fft(vector<int> &a, int invert){
	int n = a.size();
	// this is the iterative inplace implementation of FTT/NTT 
	// We need to sort the input elements according to order of recursion
	// Bit reversal permutation will be used. Refer citation
	for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1)
            j ^= bit;
        j ^= bit;

        if (i < j)
            swap(a[i], a[j]);
    }
	// Now base cases have been set in the correct order
	// Now we iterate from size of sub problem
	// Base Case: l = 1 is equal to coefficients
	// So we start from l = 2
	for (int l = 2; l<=n; l<<=1){
		// we need to set up the root of unity for n = l

		// We use the property that if w is nth root of unity then w^2 is n/2 root of unity
		// initially w_cur is log2(limit) root of unity (23rd)
		// the loop will run the number of times needed to go from 2^23 to l 
		// and each time w_cur is getting squared so at the end we get the lth root of unity
		int v = MOD-1;
		v /= l;
		int w_cur = invert ? power(primitive_i, v) : power(primitive, v);
		

		// j is the starting index
		for (int j=0; j<n; j+=l){
			int w = 1;
			for (int i=0; i<l/2; ++i){
				int first = a[i+j], second = (int)(1LL*w*a[i+j+l/2] % MOD);
				a[i+j] = first + second < MOD ? first + second : first + second - MOD;
				a[i+j+l/2] = first - second >= 0 ? first - second : first - second + MOD;

				w = (int) (1LL*w*w_cur % MOD);
			}
		}
	}

	if (invert){
		int n_i = inverse(n);

		for (int i=0; i<a.size(); ++i){
			a[i] = (int)(1LL*a[i]*n_i % MOD);
		}
	}
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
	vector<int> left(mid - l + 2);
	vector<int> right(r - (mid + 1) + 2);
	get_coefficients(left, l, mid);
	get_coefficients(right, mid+1, r);
	// Below we need to multiply left with right to get the final polynomial
	// naive implementation will be O(r^2) 
	// using FFT with NTT we will bring it down to O(rlogr)
	// using naive now
	int n = left.size(), m = right.size();
	int k = log2(max(n, m)) + 2;
	int siz = pow(2, k);

	left.resize(siz);
	right.resize(siz);
	poly.resize(siz);
	for (int i=n+1; i<siz; ++i)
		left[i] = 0;
	for (int i=m+1; i<siz; ++i)
		right[i] = 0;

	fft(left, false);
	fft(right, false);
	for (int i=0; i<siz; ++i)
		poly[i] = (int)(1LL*left[i]*right[i]%MOD);

	fft(poly, true);


}

signed main(){
	cin.tie(NULL); ios_base::sync_with_stdio(false);
	int t; cin >> t;
	while (t--){
		int n, p, r; cin >> n >> p >> r;
		int fact = 1;
		for (int i=1; i<=r; ++i)
			fact = (int)(1LL*fact*i % MOD);
		fact = inverse(fact);

		vector<int> poly(r+1);
		get_coefficients(poly, 0, r-1);

		int ans = 0;
		for (int i=1; i<=r; ++i){
			int num = (power(p, (n+1)*i) - 1);
			if (num < 0)
				num += MOD;
			int denom = (power(p, i) - 1);
			if (denom < 0)
				denom += MOD;
			if (denom != 0){
				denom = inverse(denom);
				int term = (int)(1LL*num*denom % MOD);
				int tet = (int)(1LL*term*poly[i]%MOD);
				ans = (ans + tet);
				if (ans >= MOD)
					ans -= MOD;
			}
			else {
				// common ratio of GP is congruent to 1
				int tet = (int)(1LL*((n+1)%MOD)*poly[i]%MOD);
				ans = (ans + tet);
				if (ans >= MOD)
					ans -= MOD;
			}
		}
		ans = (int)(1LL*ans*fact%MOD);
		cout << (ans>=0?ans:ans+MOD) << endl;
	}

	return 0;
}