#include <bits/stdc++.h> 
#define MOD 998244353
#define primitive 3
#define primitive_i 332748118
using namespace std;

// Check this mainly https://cp-algorithms.com/algebra/fft.html
// Also refer to https://codeforces.com/blog/entry/48798

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



// goal: NTT implementation of fast fourier transform 
void fft(int *a, int n, bool invert){

	// this is the iterative inplace implementation of FTT/NTT 
	// We need to sort the input elements according to order of recursion
	// Bit reversal permutation will be used. Refer citation
	// int lg_n = 0;
	// for (int i = 1; i<n; i<<=1)
	// 	lg_n++;
	// for (int i=0; i<n; ++i){
	// 	// Basically we compare with respect to the reversed bit 
	// 	// for instance 001 has reverse 100. When sorting would have been done
	// 	// 100 would have been at the position of 001. The reversed index
	// 	// takes the value of the original index. We do that by swapping
	// 	// To ensure that this swap is done only once we impose the inequality
	// 	if (i < reverse(i, lg_n))
	// 		swap(a[i], a[reverse(i, lg_n)]);
	// }

	// we can make do without the reverse function and instead generate the reverse
	// on the fly
	// Suppose that j already contains the reverse of i. Then by to go to i+1, 
	//we have to increment i, and we also have to increment j, but in a "reversed" 
	// number system. Adding one in the conventional binary system is equivalent to flip all tailing ones into zeros and flipping the 
	//zero right before them into a one. Equivalently in the "reversed" number system,
	// we flip all leading ones, and the also the next zero. Thus j effectively
	// generates the reversed number

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

		// we need to set up the root of unity for n = l

		// Observer that MOD-1 = c*2^23, r^c is the 2^23rd root of unity
		// (refer the other link), r^(c*2^(23)/l) will be the lth root of unity as
		// as (r^c)^(2^23) = 1 = r^(c*2^(23)/l)^l. r is the primitive root 
		// or it's inverse as needed

		int v = MOD - 1; // of the form c*2^23
		v /= l;
		int w_cur = invert ? power(primitive_i, v) : power(primitive, v);
		

		// j is the starting index
		// as elements are already in recursion order
		// all we need to do is recursively fill in values in place
		// and apply the butterfly transform
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

	// dividing by n as inverse fourier transform has a division by n 

	if (invert){
		int n_i = inverse(n);

		for (int i=0; i<n; ++i){
			a[i] = (int)(1LL*a[i]*n_i % MOD);
		}
	}
}


signed main(){
	int n, m; cin >> n >> m;
	
	int k = log2(max(n, m))+2; // pad upto 2^k (here k is sufficient such that )
	// the product will be contained with 2^k terms, that is n+m <2^k
	int siz = 1 << k;
	//poly has siz terms, degree is siz - 1

	int *pola = new int[siz];
	for (int i=0; i<=n; ++i) cin >> pola[i];
	for (int i=n+1; i<siz; ++i) pola[i] = 0;


	int *polb = new int[siz];
	for (int i=0; i<=m; ++i) cin >> polb[i];
	for (int i=m+1; i<siz; ++i) polb[i] = 0;

	fft(pola, siz, false);
	fft(polb, siz, false);
	int *final = new int[siz];
	for (int i=0; i<siz; ++i){
		final[i] = (int)(1LL*pola[i]*polb[i] % MOD);
	}

	fft(final, siz, true);

	for (int i=0; i<siz; ++i){
		cout << final[i] << " ";
	}
	cout << endl;

}