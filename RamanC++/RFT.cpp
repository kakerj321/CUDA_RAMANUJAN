#include "RFT.h"
#include <windows.h>

using namespace std;

void phitab(vector<int>& phi, int N, int i)
{
        int j = i;
        phi[j] = i;
        int p = 2;
        while (p * p <= i)
        {

            if (i % p == 0)
            {

                while (i % p == 0)
                    i = i / p;
                phi[j] -= phi[j] / p;
            }
            p++;
        }
        if (i > 1)
            phi[j] -= phi[j] / i;
    
}

int gcd(int a, int b) {
   
    // GCD(0, b) == b; GCD(a, 0) == a,
    // GCD(0, 0) == 0
    if (a == 0) {
        return b;
    }
    if (b == 0) {
        return a;
    }
    // Finding K, where K is the
    // greatest power of 2 that
    // divides both a and b.
    int k = 0;
    while (((a | b) & 1) == 0) {
        a = a >> 1;
        b = b >> 1;
        k = k + 1;
    }
   
    // Dividing a by 2 until a becomes odd
    while ((a & 1) == 0) {
        a = a >> 1;
    }
    // From here on, 'a' is always odd.
    while (b != 0) {
        // If b is even, remove all
        // factor of 2 in b
        while ((b & 1) == 0) {
            b = b >> 1;
        }

        // Now a and b are both odd.Swap if
        // necessary so a <= b, then set
        // b = b - a(which is even).
        if (a > b) {
            // Swap uand v.
            int temp = a;
            a = b;
            b = temp;
        }
        b = (b - a);
    }
    
    // restore common factors of 2
    return (a << k);


}




int cqn(vector<int>& phi_t, int q, int n, vector<int>& mu)
{


    int d = q / gcd(q, n);
    return mu[d] * phi_t[q] / phi_t[d];

}



void computeMu(int max, vector<int>& mu)
{
   

    int sqr = 0;
    sqr = sqrt(max);

    for (int i = 2; i <= sqr; i++)
    {
        if (mu[i] == 1)
        {
            for (int j = i; j <= max; j += i)
                mu[j] *= -i;
            for (int j = i * i; j <= max; j += i * i)
                mu[j] = 0;
        }
    }
    for (int i = 2; i <= max; i++)
    {
        if (mu[i] == i)
            mu[i] = 1;
        else if (mu[i] == -i)
            mu[i] = -1;
        else if (mu[i] < 0)
            mu[i] = 1;
        else if (mu[i] > 0)
            mu[i] = -1;
    }
    
}



int phi(int n)
{
   
    // Initialize result as n
    int result = n;
    // Consider all prime factors
    // of nand subtract their
    // multiples from result
    int p = 2;
    while (p * p <= n)
    {
        // Check if p is a
        // prime factor.
        if (n % p == 0)
        {
            //If yes, then
            // update nand result
            while (n % p == 0)
                n = n / p;
            result -= result / p;
        }


        p++;
    }
    // If n has a prime factor
    // greater than sqrt(n)
    // (There can be at-most
    // one such prime factor)
    if (n > 1)
        result -= result / n;
    return result;

}



int yqn(vector<int>& x, int q, int n) {
 

    int N = x.size();// (sizeof(x) / sizeof(*x));
    int sum = 0;
    int l = int(N / q); // N should be a multiple of q
    for (int j = 0; j < l; j++)
    {
        sum += (x[n + j * q - 1]);
    }
    return sum;

}

#pragma optimize( "", off )
std::chrono::steady_clock::time_point pomiar()
{
    return std::chrono::high_resolution_clock::now();
}
#pragma optimize( "", on )
float rq(vector<int>& phi_t, vector<int>& x, int q, vector<vector<int> >& cqn_cache,long long & rq_part) {
   
    auto begin_rq_part = pomiar();
    int N = x.size();
    
    volatile int sum = 0;
   
  


   
    for (int n = 1; n < (q + 1); n++) {
        sum += (yqn(x, q, n) * cqn_cache[q][n]);   
        
    }

     auto end_rq_part = pomiar();

    rq_part = std::chrono::duration_cast<std::chrono::nanoseconds>(end_rq_part - begin_rq_part).count();
    float sumf = sum;
    sumf = sumf / (float)N / (float)phi_t[q];
        return sumf;

   
}

float rf(vector<int>& phi_t, vector<int>& x, int q, vector<vector<int> >& cqn_cache, long long & rf_part) {
 

    int N = x.size();
    int sum1 = 0;
    int sum2 = 0;
    int sum3 = 0;
    auto begin_rf_part = std::chrono::high_resolution_clock::now();
    for (int n = 1; n < (q + 1); n++) {
        if (cqn_cache[q][n] == 0)
        {
            // do nothing
        }
        else if (cqn_cache[q][n] == 1)
        {
            sum1 += yqn(x, q, n);
        }
        else if (cqn_cache[q][n] == -1)
        {
            sum1 += -yqn(x, q, n);
        }
        else if (cqn_cache[q][n] == phi_t[q])
        {
            sum2 += yqn(x, q, n);
        }
        else if (cqn_cache[q][n] == -phi_t[q])
        {
            sum2 += -yqn(x, q, n);
        }
        else
        {
            sum3 += yqn(x, q, n) * cqn_cache[q][n];
        }

    }
    auto end_rf_part = std::chrono::high_resolution_clock::now();
    rf_part = std::chrono::duration_cast<std::chrono::nanoseconds>(end_rf_part - begin_rf_part).count();
 

    float sumrf1 = sum1;
    float sumrf2 = sum2;
    float sumrf3 = sum3;

    float sumrf = (sumrf1 + sumrf3) / (float)phi_t[q] + sumrf2;
    sumrf = sumrf / (float)N;
    return sumrf;

   
}

void exp(int N, int K, long long & etime_all, long long & etime_computeMu, long long & etime_cqn, long long & etime_rf, long long & etime_rq, long long & etime_rq_part, long long & etime_rf_part,long long & etime_phi)
{
    auto begin = std::chrono::high_resolution_clock::now();
    long long rf_part, rq_part;

    vector<int> x(N);
    vector<int> phi_t(K);
    auto begin_phi = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < K; i++)
    {
        phitab(phi_t, K, i);
    }
    auto end_phi = std::chrono::high_resolution_clock::now();
    etime_phi = std::chrono::duration_cast<std::chrono::nanoseconds>(end_phi - begin_phi).count();
    for (int i = 0; i < x.size(); i++) { x[i] = i % 10; }
    
    std::vector < int > mu(K, 1);
    


    auto begin_computeMu = std::chrono::high_resolution_clock::now();
    computeMu(N, mu);
    auto end_computeMu = std::chrono::high_resolution_clock::now();
    etime_computeMu = std::chrono::duration_cast<std::chrono::nanoseconds>(end_computeMu - begin_computeMu).count();
   

   
    vector<vector<int> > cqn_cache(K);
    auto begin_cqn = std::chrono::high_resolution_clock::now();
    for (int i = 0; i <= N; i++) {
        cqn_cache[i].resize(K, 1);

    }
    for (int q = 1; q < x.size() + 1; q++) {


        for (int n = 1; n < (q + 1); n++) {
            
            cqn_cache[q][n] = cqn(phi_t,q, n, mu);
            
        } 
    }
    
   
    auto end_cqn = std::chrono::high_resolution_clock::now();
    etime_cqn = std::chrono::duration_cast<std::chrono::nanoseconds>(end_cqn - begin_cqn).count();


    auto begin_rq = std::chrono::high_resolution_clock::now();
    for (int q = 1; q < x.size() + 1; q++)
    {
        
        rq(phi_t, x, q, cqn_cache,rq_part);
        etime_rq_part += rq_part;
        
    }
    
    auto end_rq = std::chrono::high_resolution_clock::now();
    etime_rq = std::chrono::duration_cast<std::chrono::nanoseconds>(end_rq - begin_rq).count();

    auto begin_rf = std::chrono::high_resolution_clock::now();
    for (int q = 1; q < x.size() + 1; q++) {

        rf(phi_t, x, q, cqn_cache,rf_part);
        etime_rf_part += rf_part;

    }
    auto end_rf = std::chrono::high_resolution_clock::now();
    etime_rf = std::chrono::duration_cast<std::chrono::nanoseconds>(end_rf - begin_rf).count();

    
    auto end = std::chrono::high_resolution_clock::now();
    etime_all = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count();

    printf("%d Time measured: %.3f seconds.\n",N, etime_all * 1e-9);

}