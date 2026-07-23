#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <climits>
#include <float.h>
#include <string>
#include <vector>
#include <algorithm>
#include <random>
#include <tr1/random>
#include <sstream> 
#include <bits/stdc++.h>
#include <sys/stat.h>
#include <sys/types.h>


#ifndef NTWK_H
#define NTWK_H


class Connettwk
{
    public:
        int N;
        unsigned int net_seed;
        int N_EZs;
        double t0;
        double speed;
        int EZ_ind;
        std::string network_name;
        std::string EZ_type="single";
        std::vector<int> EZs;
        std::vector<int> K;
        std::vector<std::vector<int>> Adj;
        std::vector<std::vector<double>> W;
        std::vector<std::vector<double>> Tau;
        std::vector<double> Exc;
        
        
        
        void load_net();
        void normal_W_Tau( double p=0.95 , double speed=60 );
        void reshape_net();
        void generate_random(int n,double p, double speed, double seed, std::string symm="symmetric");
        void set_excitabilities(double E_surround, double E_EZ);
        void save_network();
        

};




#endif