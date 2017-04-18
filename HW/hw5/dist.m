%% These are the system matrices associated with
%% a pretty fancy distillation process

A=[ -10.5234    7.4126         0         0         0;
      9.5234  -16.8957    7.4126         0         0;
           0    9.4831   16.8445    7.4126         0;
           0         0    9.4319   -16.7797    -7.4126;
           0         0         0    9.3672   16.6986];
       
B=[0.1   -0.0099;
   0     -0.0126;
   0     -0.0161;
   0     -0.0204;
   0.31  -0.0257];
       

% First, check for stability and controllability!

%% LQ Optimal control
%% Put your weights here.
Q = 1e6*eye(5);
R = eye(2);

P=are(A,B*inv(R)*B',Q);
K=inv(R)*B'*P;


%% Simulate the solution
tf=1; dt=0.001;
t=0; 
x0=[-1.4;-1.4;0.5;-0.2;-0.2]; 
x=x0;
X=[];  T=[]; 
U=[]; 
while (t<=tf);
    T=[T;t];
    X=[X,x];
    u=-K*x;
    U=[U,u];

    x=x+dt.*(A*x+B*u);
    t=t+dt;

end;

 
subplot(2,1,1); 
plot(T,X);
ylabel('x');
subplot(2,1,2); 
plot(T,U);
ylabel('u'); xlabel('t');