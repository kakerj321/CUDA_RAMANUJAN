#pragma once
#include <chrono>
#include <fstream>

#include<vector>

int gcd(int a, int b);
void computeMu(int max, std::vector<int>& mu);
int cqn(int q, int n, std::vector<int>& mu);
int phi(int n);
int yqn(std::vector<int>& x, int q, int n);
float rq(std::vector<int>& phi_t, std::vector<int>& x, int q, std::vector<std::vector<int> >& cqn_cache,long long & rq_part);
float rf(std::vector<int>& phi_t, std::vector<int>& x, int q, std::vector<std::vector<int> >& cqn_cache, long long & rf_part);
void exp(int N, int K, long long& etime_all, long long& etime_computeMu, long long& etime_cqn, long long& etime_rf, long long& etime_rq, long long& etime_rq_part, long long& etime_rf_part, long long& etime_phi);
void phitab(std::vector<int>& phi, int N, int i);
int cqn(std::vector<int>& phi_t, int q, int n, std::vector<int>& mu);