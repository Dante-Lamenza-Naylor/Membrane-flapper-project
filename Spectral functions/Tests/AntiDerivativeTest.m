addpath(genpath("C:\Users\wolfi\OneDrive\Documents\Membrane flapper project"))

trial_func = @(x)sin(x);
trial_func_antider1 = @(x)-cos(x);
trial_func_antider2 = @(x)-sin(x);
N=64;
xs = cos(pi*(0:N)/(N)); % Must be Chebyshev-Gauss-Lobatto collocation points
ys = trial_func(xs)';
ys_antider1 = trial_func_antider1(xs)';
ys_antider2 = trial_func_antider2(xs)';

% Convert to spectral space
coeffs = GetChebyCoeffs(ys);
ys_approx = Cheby2Phys(coeffs);
ys_sin_approx = 2*sum(coeffs(2:end).*sin(pi*(1:N).*(1:N)'/N)');

% Take antiderivative, enforcing boundary condition =0 at x=-1
antider1_test_coeffs = GetChebyCoeffs(trial_func_antider1(xs)');
antider1_coeffs = GetChebyAntiDerCoeffs(coeffs,trial_func_antider1(-1),N+1);
ys_antider1_approx = Cheby2Phys(antider1_coeffs(1:N+1));

% Take second antiderivative, enforcing similar boundary conditions
antider2_coeffs = GetChebyAntiDerCoeffs(antider1_coeffs,trial_func_antider2(-1),N+1);
ys_antider2_approx = Cheby2Phys(antider2_coeffs(1:N+1));

figure(1)
plot(xs,ys,'.-')
hold on
plot(xs,ys_approx,'o--')
hold off
title("y")
legend("fn","approx")

figure(2)
plot(xs,ys_antider1,'.-')
hold on
plot(xs,ys_antider1_approx,'o--')
hold off
title("ady")
legend("fn","approx")

figure(3)
plot(xs,ys_antider2,'.-')
hold on
plot(xs,ys_antider2_approx,'o--')
hold off
title("a2dy")
legend("fn","approx")