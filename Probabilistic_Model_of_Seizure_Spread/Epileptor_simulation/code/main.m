
%% Setting parameters and loading the interaction matrix
clear
clc

w=1.0;
x_0=-2.3;
tau=6667.00;
patient='P1';
EZ=64;
% To load the data unzip the connectivity zip files located in data directory
connectivity_path=strcat('../data/connectivity_',patient,'/');
W=load(strcat(connectivity_path,'weights.txt'));
W_norm=normal(W);
W=w*W_norm;

[N,~] = size(W);
%% Linear stability analysis full Jacobian
x0=x_0+zeros(N,1);
x0(EZ)=-1.6;

% initial guess
X0 = zeros(6*N,1);  
X0( 0*N+1:1*N , 1 )=-1.62;
X0( 1*N+1:2*N , 1 )=-12.2;
X0( 2*N+1:3*N , 1 )=3.1;
X0( 3*N+1:4*N , 1 )=-0.844;
X0( 4*N+1:5*N , 1 )=0;
X0( 5*N+1:6*N , 1 )=-162;

% Finding the fixed point
opt = optimset('TolFun',1e-14,'TolX',1e-14);
[Xz,fval,exitflag,output,Jhat] = fsolve(@(X)epileptor6DXdotJ(X,tau,W,x0),X0,opt);
% evaluating Jacobian, its eigenvalues and eigenvectors at the fixed point
[Xdot_in,J_in]=epileptor6DXdotJ(Xz,tau,W,x0);
[v_in,e_in]=eig(J_in);
lambda_in= max( real( eig(J_in) ) );
E_in=diag(e_in);


% Linear stability analysis: Jacobian evaluated for the surrounding while EZ in the seizure state
N=N-1;

% initial guess
x0=x_0+zeros(N,1);
X0 = zeros(6*N,1);  
X0( 0*N+1:1*N , 1 )=-1.62;
X0( 1*N+1:2*N , 1 )=-12.2;
X0( 2*N+1:3*N , 1 )=3.1;
X0( 3*N+1:4*N , 1 )=-0.844;
X0( 4*N+1:5*N , 1 )=0;
X0( 5*N+1:6*N , 1 )=-162;

% Finding the fixed point
opt = optimset('TolFun',1e-14,'TolX',1e-14);
[Xz,fval,exitflag,output,Jhat] = fsolve(@(X)epileptor6DXdotJEZ(X,tau,W,x0,EZ,0.2818),X0,opt);
% evaluating Jacobian, its eigenvalues and eigenvectors at the fixed point
[Xdot_out,J_out]=epileptor6DXdotJEZ(Xz,tau,W,x0,EZ,0.2818);
lambda_out= max( real( eig(J_out) ) );


if(lambda_in>0 && lambda_out>0)
    disp('Seizure will start in the EZ and can spread to the surrounding')
elseif(lambda_in>0 && lambda_out<0)
    disp('Seizure will start in the EZ but cannot spread to the surrounding')
else
    disp('EZ is inhibited')
end
%%

function w = normal(W)
[N,~]= size(W);
u=zeros(N*(N-1)/2,1);
k=0;
for i =1:N
    for j=i+1:N
        k=k+1;
        u(k,1)=W(i,j);
    end
    
end
u=sort(u);
uu=u(floor(length(u)*0.95));
% disp(uu)
% disp(length(u)*0.9)
w=W-diag(diag(W));
idx=w>uu;
% disp(uu)
w(idx)=uu;



w=w/(max(max(w)));


end
