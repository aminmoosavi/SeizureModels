#include "connettwk.h"

void Connettwk::load_net()
{   
    int i,j;
    double a;
    std::vector<double> v;
    std::ifstream W_file;
    std::ifstream Tau_file;
    std::ifstream EZs_file;


    // Getting EZ nodes
    EZs_file.open("./files/connectivity_"+network_name+"/EZs.txt");
    N_EZs=0;
    while(EZs_file >> a)
    {
        EZs.push_back(a-1);
        N_EZs++;
    }
    EZs_file.close();


    // Getting total number of nodes
    W_file.open("./files/connectivity_"+network_name+"/weights.txt");
    j=0;
    while(W_file >> a)
    {
        j++;
    }
    N=sqrt(j);
    W_file.close();
    
    
    // Getting interaction weights
    W_file.open("./files/connectivity_"+network_name+"/weights.txt");
    i=0;
    j=0;
    v.clear();
    while(W_file >> a)
    {
        j++;
        v.push_back(a);
        if (j==N)
        {
            j=0;
            W.push_back(v);
            v.clear();
        }    
    }
    W_file.close();
    


    // Getting tract_lengths
    Tau_file.open("./files/connectivity_"+network_name+"/tract_lengths.txt");
    i=0;
    j=0;
    t0=0;
    v.clear();
    while(Tau_file >> a)
    {
        if(a>t0)
        {
            t0=a;
        }
        j++;
        v.push_back(a);
        if (j==N)
        {
            j=0;
            Tau.push_back(v);
            v.clear();
        }      
    }
    Tau_file.close();
  


}

void Connettwk::normal_W_Tau(double p,double speedd)
{
    speed=speedd;
    int i,j;
    int M=p*N*(N-1)/2;
    std::vector<double> U;
    double v,u;
    t0=t0/speed;
    for( i=0 ; i<N-1 ; i++ )
    { 
        for( j=i+1 ; j<N ; j++ )
        {
            v=W[i][j];
            U.push_back( v );
        }
    }
    sort( U.begin() , U.end() );
    u=U[M];
    
    for( i=0 ; i<N ; i++ )
    {
        W[i][i]=0;
        for( j=0 ; j<N ; j++ )
        {
            v=W[i][j];
            Tau[i][j]=Tau[i][j]/speed;
            if(v>u)
            {
                W[i][j]=1;
            }
            else
            {
                W[i][j]=v/u;
            }
        }
        
    }
}

void Connettwk::reshape_net()
{
    int i,j,k;
    std::vector<std::vector<double>> WW;
    std::vector<std::vector<double>> TTau;
    std::vector<double> v,t;
    std::vector<int> a;
    double eps=0.00000001;

    for(i=0; i<N; i++)
    {   
        v.clear();
        t.clear();
        a.clear();
        k=0;
        for(j=0; j<N; j++)
        {
            if(W[i][j]>eps)
            {
                v.push_back(W[i][j]);
                a.push_back(j);
                t.push_back(Tau[i][j]);
                k=k+1;
            }
        }
        
        WW.push_back(v);
        TTau.push_back(t);
        Adj.push_back(a);
        K.push_back(k);
        
        

    }
    W.clear();
    Tau.clear();
    W=WW;
    Tau=TTau;
    std::cout<<"Network reshaped"<<std::endl;

}

void Connettwk::generate_random(int n, double p,double speedd, double seed, std::string symm)
{
    int i,j,k;
    double r,r1,r2,delay;
    std::vector<double> v;
    std::vector<double> t;
    std::uniform_real_distribution<double> unif(0,1);
    std::mt19937 eng;
    speed=speedd;
    // seed=(unsigned)time(NULL);
    eng.seed(seed);
    net_seed=seed;
    N=n;
    t0=0;
    std::cout<<N<<" "<<p<<std::endl;
    for(i=0; i<N; i++)
    {
        v.clear();
        t.clear();

        for(j=0; j<N; j++)
        {
            r=unif(eng);
            // std::cout<<r<<std::endl;
            
            if(r<p && i!=j)
            {
                r1=unif(eng);
                r2=unif(eng);
                delay=(0.25*r2+0.75)*200/speed;
                v.push_back( (0.2*r1+0.9) * 128.00/N );
                t.push_back(delay);
                
                if(delay>t0)
                {
                    t0=delay;
                    
                }
            }
            else
            {
                v.push_back(0.0);
                t.push_back(0.0);
            }
        }
        W.push_back(v);
        Tau.push_back(t);


    }
    
    N_EZs=N/128;
    for (i=0;i<N_EZs;i++)
    {
        EZs.push_back(i);
    }
    
    
    if(symm=="symmetric")
    {
        for (i=0;i<N;i++)
        {
            for(j=i;j<N;j++)
            {
                W[j][i]=W[i][j];
            }
        }
    }
    
    std::cout<<"Random Network Generated"<<std::endl;
  
}

void Connettwk::set_excitabilities(double E_surround, double E_EZ)
{
    int i,j;
    for(i=0;i<N;i++)
    {
        Exc.push_back(E_surround);

    }
    if(EZ_type!="single")
    {
        for(i=0;i<N_EZs;i++)
        {
            j=EZs[i];
            Exc[j]=E_EZ;
        }

    }
    else
    {
        j=EZs[EZ_ind];
        Exc[j]=E_EZ;
    }
    // std::cout << E[j-1] <<std::endl;
}

void Connettwk::save_network()
{
    int i,j;
    std::ostringstream net_seed_string, N_string;
    net_seed_string<< std::fixed;
    N_string<< std::fixed;
    net_seed_string<<net_seed;
    N_string<<N;
    std::string net_dir="./files/connectivity_"+network_name+"_N-"+N_string.str()+"-"+net_seed_string.str()+"/";
    std::string W_file_name=net_dir+"weights.txt";
    std::string Tau_file_name=net_dir+"tract_lengths.txt";
    std::string EZs_file_name=net_dir+"EZs.txt";
    std::ofstream W_out_file;
    std::ofstream Tau_out_file;
    std::ofstream EZs_out_file;
    W_out_file.open(W_file_name);
    Tau_out_file.open(Tau_file_name);
    EZs_out_file.open(EZs_file_name);
    for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
        {
            W_out_file<<W[i][j]<<" ";
            Tau_out_file<<Tau[i][j]*speed<<" ";

        }
        W_out_file<<std::endl;
        Tau_out_file<<std::endl;

    }

    for(i=0;i<N_EZs;i++)
    {
        EZs_out_file<<EZs[i]+1<<std::endl;
    }
    std::cout<<"network saved to \""+ net_dir <<"\""<<std::endl;

}

