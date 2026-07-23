clc 
clear
% Number of nodes
N=2^10;
% mu0
N0=128;
% fraction of EZ nodes
epsilon=1.0/128;
Nez=floor(epsilon*N);
% probability of connectivity in ER networks
p=0.2;
h=0.01;
x0_EZ=-1.8;
x0_c=-2.061;
Eez=h*(x0_EZ-x0_c);
% average delay time in seconds
Taubar=0.0120;
w=0.01;
E=-0.0001;
seed=floor(sum(100000*clock));
[NNe,NNs,NNr,tt,dt]=Dynamics_MF(N,N0,Nez,p,w,E,Eez,Taubar,seed);
%%
figure(1)
plot(tt,NNe/N,'LineWidth',2)
hold on
plot(tt,NNs/N,'LineWidth',2)
plot(tt,NNr/N,'LineWidth',2)
hold off
xlabel('Time [s]','interpreter','latex')
legend({'Fraction of active nodes' , 'Fraction of susceptible nodes', 'Fraction of refractory nodes'},'interpreter','latex')