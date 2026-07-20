function [Xdot,J] = epileptor6DXdotJ(X,tau,W,x0)
% 6d-epileptor 
% Returns ODE system evaluated at X
% parameters for each node
% X =[ X1(1:N),Y1(1:N),Z(1:N),X2(1:N),Y2(1:N),g(1:N)];
% W = coupling matrix for the x; M/6 by M/6
% x0 = excitability parameter (x0_i for node i; same size as number of nodes)
M=length(X);
N=M/6; % Number of nodes; X includes [x1,y1,z,x2,y2,g] 

x1=X( 0*N+1:1*N , 1 );
y1=X( 1*N+1:2*N , 1 );
z= X( 2*N+1:3*N , 1 );
x2=X( 3*N+1:4*N , 1 );
y2=X( 4*N+1:5*N , 1 );
g= X( 5*N+1:6*N , 1 );

I1=3.1;
I2=0.45;
gamma=0.01;
tau0=tau;
tau2=10;

dx1=zeros(N,1);
dy1=zeros(N,1);
dz=zeros(N,1);
dx2=zeros(N,1);
dy2=zeros(N,1);
dg=zeros(N,1);

II=eye(N,N);


x1x1=zeros(N,N); y1x1=zeros(N,N); zx1=-(1/tau0) * W; x2x1=zeros(N,N);
x1y1=II;         y1y1=-II;        zy1=zeros(N,N);    x2y1=zeros(N,N);
x1z=zeros(N,N);  y1z=zeros(N,N);  zz=-(1/tau0)*II;   x2z=-0.3*II;
x1x2=zeros(N,N); y1x2=zeros(N,N); zx2=zeros(N,N);    x2x2=zeros(N,N);
x1y2=zeros(N,N); y1y2=zeros(N,N); zy2=zeros(N,N);    x2y2=-II;
x1g=zeros(N,N);  y1g=zeros(N,N);  zg=zeros(N,N);     x2g=0.002*II;

y2x1=zeros(N,N); gx1=II;
y2y1=zeros(N,N); gy1=zeros(N,N);
y2z=zeros(N,N);  gz=zeros(N,N);
y2x2=zeros(N,N); gx2=zeros(N,N);
y2y2=(-1/tau2)*II; gy2=zeros(N,N);
y2g=zeros(N,N);  gg=-gamma*II;


for i=1:N
   
    [f1,f2,df1dx1,df1dx2,df2dx2,df1dz] = df(x1(i,1),x2(i,1),z(i,1));
    
    x1z(i,i)=-1-df1dz;
    zx1(i,i)=(1/tau0) * ( 4 + sum(W(:,i)) );
    x1x1(i,i)=-df1dx1;
    x1x2(i,i)=df1dx2;
    y2x2(i,i)=(1/tau2) * df2dx2;
    y1x1(i,i)=-10*x1(i,1);
    x2x2(i,i)= 1-3*x2(i,1)^2;
   
    
    
    dx1(i,1)=y1(i,1)-f1-z(i,1)+I1;
    dy1(i,1)=1 - 5*x1(i,1)^2 -y1(i,1);
    dz(i,1)=(1/tau0) * ( 4*(x1(i,1)-x0(i,1)) - z(i,1) - x1'*W(:,i) + sum( W(:,i) ) * x1(i,1)  );
    dx2(i,1)=-y2(i,1) + x2(i,1) - x2(i,1)^3 +I2 +0.002*g(i,1) - 0.3*(z(i,1)-3.5);
    dy2(i,1)=(1/tau2)*(-y2(i,1)+f2);
    dg(i,1)= x1(i,1) - gamma*g(i,1);   
    
end




J=cat(1,cat(2, x1x1 , x1y1 , x1z , x1x2 , x1y2 , x1g ),...
        cat(2, y1x1 , y1y1 , y1z , y1x2 , y1y2 , y1g ) ,...
        cat(2, zx1  , zy1  , zz  , zx2  , zy2  , zg) ,...
        cat(2, x2x1 , x2y1 , x2z , x2x2 , x2y2 , x2g) ,...
        cat(2, y2x1 , y2y1 , y2z , y2x2 , y2y2 , y2g) ,...
        cat(2, gx1  , gy1  , gz  , gx2  , gy2  , gg));



Xdot=cat(1, dx1,dy1,dz,dx2,dy2,dg);









end


function[f1,f2,df1dx1,df1dx2,df2dx2,df1dz] = df(x1,x2,z)



if(x1<0)
    f1=x1^3 -3*x1^2;
    df1dx1= 3 * x1^2 - 6 * x1;
    df1dx2=0;
    df1dz=0;
else
    f1=x1 * (x2 - 0.6 * (z-4)^2);
    df1dx1= x2 - 0.6 * (z-4)^2;
    df1dx2=x1;
    df1dz= -x1*1.2*(z-4);
end

if(x2>=-0.25)
    df2dx2=6;
    f2=6*(x2+0.25);
else
    f2=0;
    df2dx2=0;
end




end