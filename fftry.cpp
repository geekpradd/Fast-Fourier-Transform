#include <bits/stdc++.h>
using namespace std;

typedef complex<double> cd;
const double PI = acos(-1);

//Can replace array by vectors everywhere for extra overhead of about 14% time


void fft(cd* a, int n, bool invert){
	/*-----------
	Computes FFT if invert is false, else inverse FFT
	i.e. coeffs -> values at the n roots of unity & vice versa
	it is assumed that n is a power of 2
	------------*/

	//apply bit reversal to move from rev(i) to rev(i+1)
	for(int i = 1, j = 0; i <= n; i++){
		int bit = n >> 1; //flip leading 1's, bit represents bit to be flipped
		while(j & bit){
			j ^= bit;
			bit >>= 1;
		}
		//also flip first 0 to 1
		j ^= bit;
		//now j = rev(i)
		if(i < j){
			swap(a[i],a[j]);
		}
	}


	//now do fft iteratively for 2,4,8,..,n (n is a power of 2)
	for(int len  = 2; len <= n; len <<= 1){
		double angle = 2*PI / len;
		if(invert) angle = -angle;
		cd wcur(cos(angle),sin(angle));
		for(int i = 0; i < n; i += len){
			cd w(1,0); //start
			//for each stretch of length len, apply this level of fft
			for(int j = 0; j < len/2; j++){
				cd u = a[i+j], v = w*a[i+j+len/2];
				a[i+j] = u+v;
				a[i+j+len/2] = u-v;
				w *= wcur;
			}
		}
	}

	//if invert then divide overall by n
	if(invert){
		for(int i = 0; i < n; i++)
			a[i] /= n;
	}
}

#define MOD 998244353
#define g 3
#define g_inv 332748118
//998244353 = 1 + 7 x 17 x 2^23
/*
If p = c 2^k + 1 then g^c is a 2^k -th root of unity in the MOD modulo ring.
Here g = 3 
*/

long long int raise(long long int a, long long int b){
	if(!a) return 0;
	if(!b) return 1;
	long long int result = raise(a,b/2);
	//if(result == 0) cout << a << ' ' << b << " " << result << endl;
	result = result*result %MOD;
	if(b & 1) result = result*a % MOD;
	return result;
}

int inverse(int a){
	long long int result = raise(a,MOD-2);
	return (int)result;
}

void ntt(int *a, int n, bool invert){
	//bit reversal
	for(int i = 1, j = 0; i <= n; i++){
		int bit = n >> 1; //flip leading 1's, bit represents bit to be flipped
		while(j & bit){
			j ^= bit;
			bit >>= 1;
		}
		//also flip first 0 to 1
		j ^= bit;
		//now j = rev(i)
		if(i < j){
			swap(a[i],a[j]);
		}
	}

	for(int len = 2; len <= n; len <<= 1){
		//we need wcur to be len-th root of unity
		int val = MOD-1;
		val /= len;
		int wcur = invert? (int) raise(g_inv,val) : (int)raise(g,val);
		//int wcur = (int) raise(g,val);
		//cout << wcur << endl;
		//cout << 1LL*wcur*wcur%MOD << endl;
		//if(invert) wcur = inverse(wcur);
		//cout << wcur << endl;
		for(int i = 0; i < n; i += len){
			int w = 1;
			for(int j = 0; j < len/2; j++){
				//cout << w << ' ';
				int u = a[i+j], v = (1LL*w*a[i+j+len/2])%MOD;
				a[i+j] = (int)((u + v)%MOD);
				if(a[i+j] < 0) a[i+j] += MOD;
				a[i+j+len/2] = (u-v);
				if(a[i+j+len/2] < 0) a[i+j+len/2] += MOD;
				w = (int)((1LL*w*wcur)%MOD);
				if(w < 0) w += MOD;			
			}
			//cout << endl;
		}
	}

	if(invert){
		int n_1 = inverse(n);
		for(int i = 0; i < n; i++)
			a[i] = (1LL*a[i]*n_1)%MOD;
	}

}

int multiply(cd* a, cd* b, int n1, int n2){
	//multiplies a by b and stores the result in a;
	int n = 1;
	while(n < n1 + n2 - 1)
		n <<= 1;
	//a.resize(n); b.resize(n);
	for(int i = n1; i < n; i++)
		a[i] = 0;
	for(int i = n2; i < n; i++)
		b[i] = 0;
	fft(a,n,false);
	fft(b,n,false);
	for(int i = 0; i < n; i++){
		a[i] *= b[i];
	}
	fft(a,n,true);
	for(int i = 0; i < n; i++)
		a[i] = round(a[i].real());
	return n;
}

int multiply2(int* a, int* b, int n1, int n2){
	int n = 1;
	while(n < n1 + n2 - 1)
		n <<= 1;
	for(int i = n1; i < n; i++)
		a[i] = 0;
	for(int i = n2; i < n; i++)
		b[i] = 0;
	ntt(a,n,false);
	ntt(b,n,false);
	for(int i = 0; i < n; i++){
		a[i] = (1LL*a[i]*b[i])%MOD;
	}
	ntt(a,n,true);
}


int main(){
	int nval =  1 << 21, N =  1 << 20;
	int *a = new int[nval];
	int *b = new int[nval];
	for(int i = 0; i < N; i++)
		a[i] = b[i] = 1;
	for(int i = N; i < nval; i++)
		a[i] = b[i] = 0;

	multiply2(a,b,N,N);
	cout << a[0] << " " << a[1] << ' ' << a[2] << endl;


	return 0;
}