
function [NNe,NNs,NNr,tt,dt]=Dynamics_MF(N,N0,Nez,p,w,E,Eez,Taubar,seed)
disp('system size:')
disp(strcat('N=',num2str(N)))
disp(strcat('Nez=',num2str(Nez)))
disp(' ')

rng(seed,'twister')

NNe=[]; 
NNs=[]; 
NNr=[]; 
tt=[];
% p=0.4;
w=w*p*(N0/N);

% initial conditions
Nei=0;
Nsi=N-Nez;

Ne=Nei;
Ns=Nsi;
Nr=0;
disp('initial conditions:')
disp(strcat('Ne=',num2str(Ne)))
disp(strcat('Ns=',num2str(Ns)))
disp(strcat('Ns=',num2str(Nr)))
% Model Parameters
a=0.46;
b=0.0021;
c=1.3;
d=0.05;
Taus=1611;
Taubar=Taubar*50;

mbar=12;
if(Taubar>0)
    dt=Taubar/mbar;
else
    dt=0.01;
end

z=w*b*E*Ns;
V=w*c*(E)*Ns;
zz=0;
if(z+Eez>0)

%     T=2000000;
    M=zeros(1,2);
    M(1,1)=Nez;
    M(1,2)=1;
    NM=1;
    n=[];
    Nn=0;
    
    Ne=Nez;

    i=mbar;
    while true 
        i=i+1;
        t=(i-1)*dt;

        for j=1:NM

           Mj=M(j,1);
           jj=M(j,2);
           t_tso=dt * (i-jj);
           pj=dt*g(V,Taus,t_tso,dt,d);
           if(pj>0 && Mj>0)
               Ln=binornd(Mj,pj);
               Ne=Ne-Ln;
               M(j,1)=M(j,1)-Ln;
               
               if(Ln>0)    
                  nn=zeros(1,3);
                  nn(1,1)=Ln;
                  nn(1,2)=jj;
                  nn(1,3)=i;
                  n=cat(1,n,nn);
                  Nn=Nn+1;
                  
               end
           end
        end
        [rr,~,~]=find(M(:,1)==0);
        if(nnz(rr)>0)
            M(rr,:)=[];
            NM=NM-nnz(rr);
        end
        pste=f((z+E))*dt;
        Nss=floor(p*Ne*Ns);
        if(Nss<Ns)
            z=(w*a*dt/Taus)*zz*Ns/Nss+(w*b*E*Ns);
            pste=f((z+E))*dt;
        end
        
        if(pste>0 && Ns>0)
            Nss=floor(p*N*Ne*Ns/N);
            
            if(Nss<Ns)
                
              L=binornd(Nss,pste);
%               disp(Nss)
            else
              L=binornd(Ns,pste);
            end
%             L=binornd(Ns,pste);
            Ns=Ns-L;
            Ne=Ne+L;
            if(L>0)
                MM=zeros(1,2);
                MM(1,1)=L;
                MM(1,2)=i;
                M=cat(1,M,MM);
                NM=NM+1;
            end

        end

        Nr=N-Ns-Ne;

%         disp(strcat(num2str(i),"  ",num2str(Ne)))
        zz=0;
        if(nnz(n)>0)
            KK=2*n(:,3)-n(:,2)-i+mbar;
            [rr,~,vv]=find(KK<=0);
            n(rr,:)=[];
            KK(rr,:)=[];
            Nn=Nn-sum(vv);
            kk=h(KK);
            zz=zz+dot(kk,n(:,1));
          
        end
        if(nnz(M(:,1))>0)
            JJ=h(i-M(:,2)-mbar);
            zz=zz+dot(JJ,M(:,1));
        end
        z=(w*a*dt/Taus)*zz+(w*b*E*Ns);
        V=(w*c*(E)*Ns);

        
        tt=cat(1,tt,t);
        NNe=cat(1,NNe,Ne);
        NNs=cat(1,NNs,Ns);
        NNr=cat(1,NNr,Nr);
        if(Nn+NM==0)
            break
        end
    end
else
    Ne=0;
    Ns=N-Nez;
    Nr=0;
    t=0;
    tt=t;
    NNe=Ne;
    NNs=Ns;
    NNr=Nr;
    
    
end

tt=tt/50;
end





function y=g(V,T,t_tso,dt,d)

    
    mu=T/(1.0-V);
    sig=d*mu;
%     sig=sqrt(3)*mu/(32);
    if(t_tso<mu-sig)
        y=0;
    elseif(t_tso>=mu+sig)
        y=1/dt;
    else
        y=1.0/(sig+mu-t_tso);
        if(y>1.0/dt)
           y=1.0/dt; 
        end
    end
end


function y=f(x)
    if(x>1)
        y=1;
    elseif(x<0)
        y=0;
    else
        y=x;
        
    end
end



function y=h(x)
idx=x>0;
y=x.*idx;
end

% function z=getz(M,n,mbar,i)
% 
% kk=h(2*n(:,3)-n(:,2)-i+mbar);
% jj=h(i-M(:,2)-mbar);
% z=kk'*n(:,1)+jj'*M(:,1);
% 
% end
