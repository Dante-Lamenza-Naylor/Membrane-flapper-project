% Note to self: addpath function takes a lot of compute time...

N = 64;
Ae = 1;
R = 1;
U = 1;
sigma = 2*pi/U;
lambda_0 = 1.1;
amp = 1;
k_edge = 1;

y_targets = zeros(N+1,1);
% dydt0     = zeros(N+1,1);
dydt0 = amp*sigma*ones(N+1,1);
x = cos(pi*(0:N)/N)'; % Should never change...

params = {"N", N;
          "Ae",Ae;
          "R" , R;
          "sigma",sigma;
          "lambda_0",lambda_0;
          "amp",amp;
          "U",U;
          "k_edge",k_edge;
          "y_targets",y_targets;
          "x",x};

tmax = 2;
dt = 1e-3;
t0 = 0;

basket_a = solver_basket(dt);
basket_a.SetUpMemory(dt,dt,tmax,6);

outs = timeintzero(@(t,y)RigidWingSolver(t,y,basket_a,params),t0,dt,tmax,[y_targets ; dydt0]);

%%
RigidWingViz(basket_a,params,1,"N=64,dt=1e-5")
