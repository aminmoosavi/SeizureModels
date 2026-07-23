#include "connettwk.h"
#include "dynamics.h"
#include <iomanip>
#include <chrono>
#include <ratio>
#include <ctime>

int main(int argc, char *argv[])
{   
    const std::string network_name=argv[1];
    const std::string EZindex=argv[2];
    const std::string w_str=argv[3];
    const std::string E_str=argv[4];
    const std::string Eez_str=argv[5];
    const std::string NR_str=argv[6];
    
    int EZ_ind=std::stoi(EZindex)-1;
    double wi=std::stod(w_str);
    double Ei=std::stod(E_str);
    double Eez=std::stod(Eez_str);
    int NR=std::stoi(NR_str);
    std::cout << "wi=" << wi << " Ei=" << Ei <<" Eez=" << Eez <<" NR=" << NR <<std::endl;
    int N;
    std::vector<int> susceptibles;
    std::vector<int> exciteds;
    std::vector<double> Exc;
    unsigned long int seed;
    int i,j,k,l;
    long double x,y;
    double a=0.46,b=0.0021,c=1.3,d=0.05,H,T=1611,w,E;
    double dw=0.075,dE=0.0115;
    int Nw=1, NE=1;
    H=b/0.21;
    std::ostringstream w_string, E_string , EZ_string, Tau_string;
    w_string<< std::fixed;
    EZ_string<< std::fixed;
    E_string<< std::fixed;
    Tau_string<< std::fixed;
    w_string<< std::setprecision(6);
    E_string<< std::setprecision(6);

    std::cout << "Minimal Modle" << std::endl;
    std::cout << "Network : " << network_name <<std::endl;
    std::cout << "EZ index : " << EZ_ind <<std::endl;
    std::cout << "a="<<a<<" b="<<b<<" c="<<c<<" d="<<d<<" tau="<<T<<std::endl;
    
    

    for(i=0;i<Nw;i++)
    {
        w=wi+i*dw;
        for(j=0;j<NE;j++)
        {
            E=Ei+j*dE;
            std::cout << "w=" << w << " E=" << E <<std::endl;
            std::cout << "i=" << i << " j=" << j <<std::endl;
            std::cout << "Nw=" << Nw << " NE=" << NE <<std::endl;
            Dynamics Y;
            std::cout<<"simulation object created"<<std::endl;
            Y.set_parameters(a,b,c,d,H,T,w,E);
            Y.network_name=network_name;
            
            Y.load_net();
            std::cout<<"network loaded"<<std::endl;
            Y.normal_W_Tau(0.95,60.0);
            std::cout<<"network normalized"<<std::endl;
            Y.reshape_net();
            Y.EZ_ind=EZ_ind;
            std::cout<<"EZ node = "<<Y.EZs[Y.EZ_ind]+1<<std::endl;
            Y.set_excitabilities(E,Eez);
            std::cout<<"Set excitabilities"<<std::endl;
            Y.initialize();
            std::cout<<"simulation initialized"<<std::endl;
                    
            
            w_string.str(std::string());
            E_string.str(std::string());
            EZ_string.str(std::string());
            Tau_string.str(std::string());


            w_string<<w;
            E_string<<E;
            EZ_string<<Y.EZs[Y.EZ_ind]+1;
            Tau_string<<T;

            std::string out_dir="./data/"+network_name+"/";

            std::string file_name_order=out_dir+"spread_w-"+w_string.str()+"_E"+E_string.str()+".txt";
            std::string file_name_onsets=out_dir+"onsets_w-"+w_string.str()+"_E"+E_string.str()+".txt";
            std::string file_name_offsets=out_dir+"offsets_w-"+w_string.str()+"_E"+E_string.str()+".txt";
            std::cout<<file_name_order<<std::endl;
            std::ofstream out_file_order;
            std::ofstream out_file_onsets;
            std::ofstream out_file_offsets;

            out_file_order.open(file_name_order);
            out_file_onsets.open(file_name_onsets);
            out_file_offsets.open(file_name_offsets);


            for (k=0;k<NR;k++)
            {
                std::chrono::high_resolution_clock::time_point t1;
                t1 = std::chrono::high_resolution_clock::now();
                seed=std::chrono::duration_cast<std::chrono::nanoseconds>(t1.time_since_epoch()).count();

                // std::cout<<"seed="<<seed<<std::endl;
                
                Y.seed=seed;
                Y.run();
                out_file_order<< Y.order <<std::endl;
                // std::cout<<Y.order<<std::endl;
                for (l=0;l<Y.N;l++)
                {
                    out_file_onsets<< Y.onsets[l]/50.0<<" ";
                    out_file_offsets<< Y.offsets[l]/50.0<<" ";
                    
                }
                out_file_onsets<<std::endl;
                out_file_offsets<<std::endl;
                
            }





        }

    } 


    
    





    
    

   

    return 0;
}