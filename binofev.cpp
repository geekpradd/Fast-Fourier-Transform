#include <bits/stdc++.h>
#define int long long
#define primitive 3
#define primitive_i 332748118
#define MOD 998244353
// primitive_i is the inverse of the primitive root of MOD
using namespace std;

int power(int a, int b){
	if (b==0)
		return 1;
	int num = power(a, b/2);
	num = (num*num)%MOD;
	if (b%2)
		return num*a%MOD;
	else
		return num;
}
int inverse(int a){
	return power(a, MOD-2);
}


// goal: NTT implementation of fast fourier transform 
void fft(int *a, int n, int invert){

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

		// Observer that MOD-1 = c*2^23, r^c is the 2^23rd root of unity
		// (refer the other link), r^(c*2^(23)/l) will be the lth root of unity as
		// as (r^c)^(2^23) = 1 = r^(c*2^(23)/l)^l. r is the primitive root 
		// or it's inverse as needed

		int v = MOD-1;
		v /= l;
		int w_cur = invert ? power(primitive_i, v) : power(primitive, v);
		

		// j is the starting index
		// as elements are already in recursion order
		// all we need to do is recursively fill in values in place
		// and apply the butterfly transform
		for (int j=0; j<n; j+=l){
			int w = 1;
			for (int i=0; i<l/2; ++i){
				int first = a[i+j], second = w*a[i+j+l/2] % MOD;
				a[i+j] = first + second < MOD ? first + second : first + second - MOD;
				a[i+j+l/2] = first - second >= 0 ? first - second : first - second + MOD;

				w = w*w_cur % MOD;
			}
		}
	}

	// dividing by n as inverse fourier transform has a division by n 
	if (invert){
		int n_i = inverse(n);

		for (int i=0; i<n; ++i){
			a[i] = a[i]*n_i % MOD;
		}
	}
}

void get_coefficients(int *poly, int siz, int l, int r){
	if (l==r){
		poly[0] = -1*(l); poly[1] = 1;
		if (poly[0]!=0)
			poly[0] = poly[0] + MOD;
		return;
	}
	int mid = (l+r)/2;
	int *left = new int[siz]; 
	int *right = new int[siz];
	for (int i=0; i<siz; ++i)
		left[i] = right[i] = 0;
	int l_cap = log2(mid-l+1) + 1; int l_siz = 1 << l_cap;
	int r_cap = log2(r-mid) + 1; int r_siz = 1 << r_cap;

	get_coefficients(left, l_siz,  l, mid);
	get_coefficients(right, r_siz, mid+1, r);
	// Below we need to multiply left with right to get the final polynomial
	// naive implementation will be O(r^2) 
	// using FFT with NTT we will bring it down to O(rlogr)
	

	fft(left, siz, false);
	fft(right, siz, false);
	for (int i=0; i<siz; ++i)
		poly[i] = left[i]*right[i]%MOD;

	fft(poly, siz, true);


}

signed main(){
	cin.tie(NULL); ios_base::sync_with_stdio(false);
	int t; cin >> t;
	int *facto = new int[100001];
	facto[1] = 1;
	for (int i=2; i<=100000; ++i)
		facto[i] = facto[i-1]*i % MOD;
	while (t--){
		int n, p, r; cin >> n >> p >> r;
		
		int fact = inverse(facto[r]);

		int cap = log2(r) + 1; int siz = 1 << cap;
		
		int *poly = new int[siz];
		for (int i=0; i<siz; ++i) poly[i] = 0;
		get_coefficients(poly, siz, 0, r-1);

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
				int term = num*denom % MOD;
				int tet = term*poly[i]%MOD;
				ans = (ans + tet);
				if (ans >= MOD)
					ans -= MOD;
			}
			else {
				// common ratio of GP is congruent to 1
				int tet = ((n+1)%MOD)*poly[i]%MOD;
				ans = (ans + tet);
				if (ans >= MOD)
					ans -= MOD;
			}
		}
		ans = ans*fact%MOD;
		cout << (ans>=0?ans:ans+MOD) << endl;
	}

	return 0;
}