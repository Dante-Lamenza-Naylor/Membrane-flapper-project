function p = AnalyticWingPressure(xN,t,beta0,beta1,sigma,U)
% From Moore 2014
% Calculates dimensionless pressure across membrane wing for small
% amplitude heaving and pitching motions
% xN: amount of nodes in membrane domain
% t: current time
% beta0 and beta1: movement parameters
% Take sigma = 2pi/U

x = cos(pi*(0:xN)/xN)';

% Reduced frequency
thisTheo = Theodorsen(sigma);

% Pressure paramters
a0 = 0* (-2*pi*1i*U*thisTheo*beta0...
    + 2*pi*1i*U*(1-thisTheo)*beta1...
    - 2*U^2*thisTheo*beta1)/2;

a1 = (2*pi^2*beta0 - 4*pi*1i*U*beta1);

a2 = pi^2*beta1;

% Dimensionless pressure drop across wing
p = real((a0*sqrt((1-x)./(1+x)) + 2*a1*sqrt(1-x.^2) + 4*a2*x.*sqrt(1-x.^2)).*exp((2*pi*t - pi/2)*1i));