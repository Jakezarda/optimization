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


int main() {
    
    std::random_device seeder;
    std::mt19937_64 gen(seeder());
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    double start;
    double end;
    int N_trials = 1000;
    int N_darts = 1048576;
    std::vector<double> pis;
    
    start = omp_get_wtime();          
                   
    for (int i = 0; i < N_trials; ++i) {
        int count = 0;  
        for (int j = 0; j < N_darts; ++j){
            double x = dist(gen);
            double y = dist(gen);
            double r = std::sqrt(x*x + y*y);
        
            if (r <= 1.0) {                 // ask how to get c++ to output to a file
                count++;
            }
        }
        double pi = 4.0*int(count)/int(N_darts);
        pis.push_back(pi);
    }
    
    end = omp_get_wtime();
    
    double average;
    double sumtotal = 0;
    for (int i = 0; i < pis.size(); ++i) {
        sumtotal += pis[i];
    }
    
    average = sumtotal/pis.size();
    double variance = 0;
    double stdev = 0;
    
    for (int i = 0; i < pis.size(); ++i){
        variance += pow(pis[i] - average,2);
    }
    
    variance = variance/pis.size();
    stdev = sqrt(variance);
        
    std::cout.precision(6);
    std::cout << "Time to compute " << N_trials << " estimates of pi: " <<  end - start << "s" << std::endl;
    std::cout << "Average found for pi was: " << average << std::endl;
    std::cout << "Standard deviation is: " << stdev << std::endl;
    std::ofstream fout("PiCPP.dat");
    fout.precision(15);
    for(int i = 0; i < N_trials; ++i){
        fout << pis[i] << "\n";   
    }
    
    
    return 0;
}






