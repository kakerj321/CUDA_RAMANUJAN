#include <iostream>
#include <algorithm>
#include <stdio.h>
#include <math.h>
#include <cstdlib>
#include <vector>
#include "RFT.h"
#include <chrono>
#include <fstream>


using namespace std;



int main()
{
    std::fstream time_all, time_computeMu, time_cqn, time_rf, time_rq, time_rf_part, time_rq_part,time_phi;
    time_all.open("time_all.csv", std::ios::app);
    time_computeMu.open("time_computeMu.csv", std::ios::app);
    time_cqn.open("time_cqn.csv", std::ios::app);
    time_rf.open("time_rf.csv", std::ios::app);
    time_rq.open("time_rq.csv", std::ios::app);
    time_rf_part.open("time_rf_part.csv", std::ios::app);
    time_rq_part.open("time_rq_part.csv", std::ios::app);
    time_phi.open("time_phi.csv", std::ios::app);

    
    time_all << "Numbers of elements" << ";" << "Time measured" << "\n";
    time_computeMu << "Numbers of elements" << ";" << "Time measured" << "\n";
    time_cqn << "Numbers of elements" << ";" << "Time measured" << "\n";
    time_rf << "Numbers of elements" << ";" << "Time measured" << "\n";
    time_rq << "Numbers of elements" << ";" << "Time measured" << "\n";
    time_rf_part << "Numbers of elements" << ";" << "Time measured" << "\n";
    time_rq_part << "Numbers of elements" << ";" << "Time measured" << "\n";
    time_phi << "Numbers of elements" << ";" << "Time measured" << "\n";

    for (int N = 100; N <= 20000; N += 100)
    {
        int K = N + 1;
        long long etime_all=0, etime_computeMu=0, etime_cqn=0, etime_rf=0, etime_rq=0, etime_rq_part=0, etime_rf_part=0,etime_phi=0;
        exp(N, K, etime_all, etime_computeMu,  etime_cqn, etime_rf, etime_rq, etime_rq_part,etime_rf_part,etime_phi);

        time_all << N << ";" << etime_all * 1e-9 << "\n";
        time_computeMu << N << ";" << etime_computeMu * 1e-9 << "\n";
        time_cqn << N << ";" << etime_cqn * 1e-9 << "\n";
        time_rf << N << ";" << etime_rf * 1e-9 << "\n";
        time_rq << N << ";" << etime_rq * 1e-9 << "\n";
        time_rf_part << N << ";" << etime_rf_part * 1e-9 << "\n";
        time_rq_part << N << ";" << etime_rq_part * 1e-9 << "\n";
        time_phi << N << ";" << etime_phi * 1e-9 << "\n";
    }
   
    time_all.close();
    time_computeMu.close();
    time_cqn.close();
    time_rf.close();
    time_rq.close();
    time_rf_part.close();
    time_rq_part.close();
    time_phi.close();

    return EXIT_SUCCESS;
}
