% Note to self: addpath function takes a lot of compute time...

N = 64;
Ae = 1;
R = 1;
U = 1;
sigma = 2*pi/U;
lambda_0 = 1.1;
amp = 1;
k_edge = 1;

x = flip(cos(pi*(0:N)/N)'); % Should never change...
y_targets = -amp*cos(x*pi/2);
dydt0     = zeros(N+1,1);

%% Linear analytical limit
H_m = 2;
a = 4*R;
b = 2*Ae*(lambda_0 - 1);
natfreq_lin = sqrt(b/a)/(2*H_m);
natper_lin  = 1/natfreq_lin;

%% Nonlinear limit
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

tmax = natper_lin/2;
dt = 1e-3;
t0 = 0;

basket_a = solver_basket(dt);
basket_a.SetUpMemory(1e-4,1e-3,tmax,9);


outs = timeintzero(@(t,y)StillVacuumSolver(t,y,basket_a,params),t0,dt,tmax,[y_targets ; dydt0]);

% Analytical comparison
t_sol           = cell2mat(basket_a.LTM(1,:));
y_sol           = cell2mat(basket_a.LTM(2,:));
lin_sol = amp * cos(sqrt(b/a) * pi * t_sol(1:end-1)/H_m).*sin(pi*(x-1)/2);

error = max(abs(y_sol - lin_sol),[],'all');
%%
loglog(amps,errors,'o')
xlabel("Initial amplitude $a$",'interpreter','latex')
ylabel("Max error",'interpreter','latex')
%%
StillVacuumViz(basket_a,params,1,"test3")
