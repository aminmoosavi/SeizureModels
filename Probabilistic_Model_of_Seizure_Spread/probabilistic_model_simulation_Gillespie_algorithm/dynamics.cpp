#include "dynamics.h"
#include "connettwk.h"



double Dynamics::f(double x)
{
    if(x>1)
    {
        return 1;
    }
    else if(x<0)
    {
        return 0;
    }
    else
    {
        return x;
    }
}

// double Dynamics::f(double x)
// {
//     if(x<0)
//     {
//         return 0;
//     }
//     else
//     {
//         return x;
//     }
// }

double Dynamics::g(double t_tso , double yd)
{
    double mu, sig;
    mu=T/(1.0-yd);
    // sig=sqrt(3)*mu/(32);
    sig=d*mu;
    if(t_tso<mu-sig)
    {
        return 0;
    }
    else if(t_tso>=mu+sig)
    {
        return 1/h;
    }
    else
    {
        if(1.0/(sig+mu-t_tso)>1/h)
        {
            return 1/h;
        }
        else
        {
            return 1.0/(sig+mu-t_tso);
        }
        
    }
}

void Dynamics::set_parameters
(double a, double b, double c,double d, double H, double T,double w,double E)
{
    this->a=a;
    this->b=b;
    this->c=c;
    this->d=d;
    this->H=H;
    this->T=T;
    this->w=w;
    this->E=E;
}

void Dynamics::initialize()
{
    int i,j,k,l;
    int EZ,EZZ;
    for(i=0; i<N ;i++)
    {

        status.push_back(-1);
        susceptibles.push_back(i);
        exciteds.push_back(N);
        refractories.push_back(N);
        P.push_back(0);
        y.push_back(0);
        yD.push_back(0);
        Delta.push_back(0);
        DeltaD.push_back(0);
        onsets.push_back(0);
        offsets.push_back(0);
        durations.push_back(0);
        n_sus=N;
        n_exc=0;
        n_ref=0;
        t=t0;
        t_seiz_off_last=0;
        t_seiz_on_last=0;
        t_transition=0;
    } 
                       
}


void Dynamics::EZ_activation()
{
    int i,j,l=0;
    double lambda,y1;



    if(EZ_type=="single")
    {
        j=EZs[EZ_ind];
        y1=y[j];
        y[j]=this->history(j);
        Delta[j]=y[j]-y1;
        lambda=this->f(y[j]);
        
        if(y[j]>0.00)
        {
            status[j]=1;
            exciteds[0]=j;
            onsets[j]=t;
            t_seiz_on_last=t;
            t_transition=t;
            n_exc++;
            susceptibles[j]=susceptibles[n_sus-1];
            susceptibles[n_sus-1]=N;
            n_sus--;
            
        }
    }
    else if(EZ_type!="single")
    {
        for(i=0;i<N_EZs;i++)
        {
            j=EZs[i];
            y1=y[j];
            y[j]=this->history(j);
            Delta[j]=y[j]-y1;
            lambda=this->f(y[j]);
            if(lambda>0)
            {
                // std::cout<<"multiple EZ nodes"<<std::endl;
                l=1;
                break;
            }
        }
        if(l==1)
        {

            for(i=0;i<N_EZs;i++)
            {
                j=EZs[i];
                status[j]=1;
                exciteds[i]=j;
                onsets[j]=t;
                t_seiz_on_last=t;
                t_transition=t;
                n_exc++;
                susceptibles[j]=susceptibles[n_sus-1];
                susceptibles[n_sus-1]=N;
                n_sus--; 
            }

        }
    }
}


void Dynamics::one_step()
{   
    int i,j,l=0;
    double tt,r,S1=0,lambda,y1,y2,endtime,t11;
    
    std::uniform_real_distribution<double> unif(0,1);
    std::exponential_distribution<double> exponential_dist(1.0);
    
    t=t+h;
    t11=t;
    for(i=0;i<n_ref;i++)
    {
        j=refractories[i];
        if(status[j]==0 && t>offsets[j]+durations[j])
        {
            status[j]=-2;
            t_transition=offsets[j]+durations[j];
            if(t_transition<t11)
            {
                t11=t_transition;
            }
            t_transition=t11;

        }
    
    }
    t=t11;
    


    for(i=0;i<n_exc;i++)
    {
        j=exciteds[i];
        tt=(t-onsets[j]);
        // std::cout<<j<<" "<<t<<std::endl;
        
        
        if(t>t_transition+t0+2*h)
        {
            y[j]=y[j]+Delta[j];
            yD[j]=yD[j]+DeltaD[j];
        }
        else
        {
            y1=y[j];
            y[j]=this->history(j);
            Delta[j]=y[j]-y1;
            y2=yD[j];
            yD[j]=this->historyD(j);
            DeltaD[j]=yD[j]-y2;
        }

        
        lambda=this->g( tt , yD[j] );
        S=S+lambda*h;
        S1=S1+lambda*h;
        P[l]=lambda*h;
        l=l+1;
    }
    for(i=0;i<n_sus;i++)
    {
        j=susceptibles[i];
        if(t>t_transition+t0+2*h)
        {
            y[j]=y[j]+Delta[j];
            yD[j]=yD[j]+DeltaD[j];
        }
        else
        {
            y1=y[j];
            y[j]=this->history(j);
            Delta[j]=y[j]-y1;
            y2=yD[j];
            yD[j]=this->historyD(j);
            DeltaD[j]=yD[j]-y2;
            
        }
        lambda=this->f(y[j]);
        S=S+lambda*h;
        S1=S1+lambda*h;
        P[l]=lambda*h;
        l=l+1;
    }
    
    if(S>t_silent)
    {
        t_transition=t+(t_silent-S)/S1;
        r=unif(eng);
        i=this->random_select(l,r,S1);
        if(i<n_exc)
        {
            j=exciteds[i];
            exciteds[i]=exciteds[n_exc-1];
            exciteds[n_exc-1]=N;
            n_exc--;
            refractories[n_ref]=j;
            n_ref++;
            offsets[j]=t_transition;
            status[j]=0;
            durations[j]=offsets[j]-onsets[j];
            t_seiz_off_last=t_transition;
            duration_last=offsets[j]-onsets[j];
        }
        else
        {
            i=i-n_exc;
            j=susceptibles[i];
            susceptibles[i]=susceptibles[n_sus-1];
            susceptibles[n_sus-1]=N;
            n_sus--;
            exciteds[n_exc]=j;
            n_exc++;
            onsets[j]=t_transition;
            status[j]=1;
            t_seiz_on_last=t_transition;


        }
        t=t_transition;
        t_silent=exponential_dist(eng);
        S=0;
       


    }
       
}

void Dynamics::run()
{   
    int i;
    std::exponential_distribution<double> exponential_dist(1.0);
    
    
    
    eng.seed(seed);

    n_sus=N;
    n_exc=0;
    n_ref=0;
    t=t0;
    t_seiz_off_last=0;
    t_seiz_on_last=0;
    duration_last=0;
    t_transition=0;
    for(i=0;i<N;i++)
    {
        status[i]=-1;
        susceptibles[i]=i;
        exciteds[i]=N;
        refractories[i]=N;
        onsets[i]=0;
        offsets[i]=0;
        durations[i]=0;
    }
    
    this->EZ_activation();
    S=0;
    if(n_exc>0)
    {
        t_silent=exponential_dist(eng);
        // std::cout<<"t_silent= "<<t_silent<<std::endl;
        

        
        while(true)
        {
            
            this->one_step();

            if(n_exc==0 && t>t0+t_seiz_off_last+duration_last)
            {
                
                break;
            }
            

        }

    }
    

    order=N-n_sus-n_exc;
    // std::cout<<"order="<<order<<std::endl;
    // std::cout<<offsets[EZs[EZ_ind]]-onsets[EZs[EZ_ind]]<<std::endl;


}


double Dynamics::history(int n)
{
    int i,j;
    double tt,ttt,y1,duration;
    
    y1=Exc[n];
    for(i=0; i<K[n];i++)
    {
        j=Adj[n][i];
        if(status[j]==-1)
        {
            y1=y1+b*w*Exc[j]*W[n][i];
        }
        else if(status[j]==1)
        {
            tt=t-onsets[j]-Tau[n][i];
            if(tt>0)
            {
                y1=y1+a*(tt/T  )*w*W[n][i];
            }
            else if(tt<0)
            {
                y1=y1+b*w*Exc[j]*W[n][i];
            }
        }
        else if(status[j]==0)
        {
            ttt=t-offsets[j]-Tau[n][i];
            tt=t-onsets[j]-Tau[n][i];
            duration=offsets[j]-onsets[j];
            if(ttt<0)
            {
                y1=y1+a*(tt/T )*w*W[n][i];
            }
            else if(ttt>0 && ttt<duration)
            {
                y1=y1+a*( (duration-ttt)/(T) )*w*W[n][i];
            }

        }

    }
    return y1;
}





double Dynamics::historyD(int n)
{
    int i,j;
    double tt,y1;
    double x0b;
    x0b=0;
    
    y1=0;
    for(i=0; i<K[n];i++)
    {
        j=Adj[n][i];
        if(status[j]==-1)
        {
            y1=y1+c*w*( Exc[j]+x0b )*W[n][i];
        }
        else if(status[j]==1)
        {
            tt=t-onsets[j]-Tau[n][i];
            if(tt<0)
            {
                y1=y1+c*w*( Exc[j]+x0b )*W[n][i];
            }
            // else if(tt>0 && tt<T)
            // {
            //     y1=y1+c*w*( Exc[j]+x0b )*W[n][i]*(1-tt/T);
            // }
        }
       
    }
    return y1;
}



int Dynamics::random_select(int l,double r,double S1)
{
    int i,j;
    double Cdf=0;
    for(i=0;i<l;i++)
    {
        Cdf=Cdf+P[i]/S1;
        // std::cout<<Cdf<<std::endl;
        if(Cdf>r)
        { 
            break;
        }
    }
    
    return i;
}
