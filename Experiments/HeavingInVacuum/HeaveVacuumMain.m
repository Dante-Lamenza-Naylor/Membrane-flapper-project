% Note to self: addpath function takes a lot of compute time...

N = 64;
Ae = 1;
R = 1;
lambda_0 = 1.1;

nat_freq_est = (pi/sqrt(8))*sqrt(Ae*(lambda_0-1)/(R+0.5)); % From Gali

sigma = 1*nat_freq_est;
U = 2*pi/sigma;
amp = 1;
k_edge = 1;

x = flip(cos(pi*(0:N)/N)'); % Should never change...
y_targets = zeros(N+1,1);
dydt0     = zeros(N+1,1);

params = {"N", N;                % 1
          "Ae",Ae;               % 2
          "R" , R;               % 3
          "sigma",sigma;         % 4
          "lambda_0",lambda_0;   % 5
          "amp",amp;             % 6
          "U",U;                 % 7
          "k_edge",k_edge;       % 8
          "y_targets",y_targets; % 9
          "x",x};                % 10



tmax = 20;
dt = 1e-3;
t0 = 0;

basket_a = solver_basket(dt);
basket_a.SetUpMemory(dt,dt,tmax,8);

outs = timeintzero(@(t,y)HeaveVacuumSolver(t,y,basket_a,params),t0,dt,tmax,[y_targets ; dydt0]);

%%
HeaveVacuumViz(basket_a,params,1,"Sigma17")
