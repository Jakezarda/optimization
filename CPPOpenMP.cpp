#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <random>
#include <omp.h>

struct parameters{
    int N, threads;
};

parameters parser(int argc, char *argv[]) {
    parameters p;
    p.N = 1000;                                     //at least 1000 points
    p.threads = omp_get_max_threads();
    for (int i = 1; i < argc; ++i) {
        std::string temp(argv[i]);
        if (temp == "-N") {
            p.N = std::stoi(argv[i + 1]);
            i++;
        } else if (temp == "-threads") {
            p.threads = std::stoi(argv[i + 1]);
        }
    }
    return p;
}

int main(int argc, char *argv[]) {
    parameters p = parser(argc, argv);
    int N_darts = 1024*1024
    int N_trials = 1000
    //rng
    //NO CLUE WHAT THE PARAMETERS THING IS
    std::vector<double> pis(p.threads);
    
#pragma omp parallel
    {
        int tid = omp_get_thread_num();
        std::random_device seeder;
        std::mt19937_64 gen(seeder());
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        
        std::vector <int> count(omp_get_max_threads());
        for (int i = 0; i < p.N; ++i) {
            int tid = omp_get_thread_num();          //gets thread ID
            double x = dist(gen);
            double y = dist(gen);
            double r = std::sqrt(x*x + y*y);
            
            if (r <= 1.0) {
                count[tid]++;
            }
        }
        
        pis[] = 4.0*double(count)/double(p.N);
    }
    
    double pi = 0.0;

    for (int i = 0; i < p.threads; ++i){
        pi += pis[i]   
        std::cout << "pis[" << i << "] = " <<pis[i] << std::endl;
    }
    
    pi /= double(p.threads);
    
    
    std::cout.precision(15);
    std::cout << "pi_est = " << pi << std::endl;
    std::cout << "pi_real = " << M_PI << std::endl;
    double percent_error = 100.0*std::abs(pi - M_PI)/M_PI;
    std::cout.precision(6);
    std::cout << "Percent Error: " << percent_error << "%\n";
    
    return 0;
}
    
