#include "connettwk.h"


#ifndef DYNAMICS_H
#define DYNAMICS_H

class Dynamics : public Connettwk
{   
    double t,h=0.05,S,t_seiz_off_last,t_seiz_on_last,t_transition,duration_last;  
    std::vector<int> susceptibles;
    std::vector<int> exciteds;
    std::vector<int> refractories;
    std::vector<double> P;
    std::vector<double> y;
    std::vector<double> yD;
    std::vector<double> Delta;
    std::vector<double> DeltaD;
    int n_sus;
    int n_exc;
    int n_ref;
    double f(double x);
    double g(double t_ts0 , double yd);
    int random_select(int l,double r,double S1);
    double t_silent;
    std::mt19937_64 eng;
    
    

    public:
        double a;
        double b;
        double c;
        double d;
        double H;
        double T;
        double w;
        double E;
        unsigned int seed;
        std::vector<int> status;
        std::vector<double> onsets;
        std::vector<double> offsets;
        std::vector<double> durations;
        unsigned int order; 
        
        
        void set_parameters(double,double,double,double,double,double,double,double);
        void initialize();
        void one_step();
        void run();
        void EZ_activation();
        double history(int i);
        double historyD(int i);
        // double history(int i);
                
};


#endif