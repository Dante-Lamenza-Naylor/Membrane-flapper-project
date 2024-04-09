trial_func = @(x)sin(x);
trial_func_der1 = @(x)cos(x);
trial_func_der2 = @(x)-sin(x);
N=64;
xs = cos(pi*(0:N)/(N)); % Must be Chebyshev-Gauss-Lobatto collocation points
ys = trial_func(xs)';
ys_der1 = trial_func_der1(xs)';
ys_der2 = trial_func_der2(xs)';

% Convert to spectral space
coeffs = GetChebyCoeffs(ys);
ys_approx = Cheby2Phys(coeffs);
ys_sin_approx = 2*sum(coeffs(2:end).*sin(pi*(1:N).*(1:N)'/N)');

% Take antiderivative, enforcing boundary condition =0 at x=-1
der1_test_coeffs = GetChebyCoeffs(trial_func_der1(xs)');
der1_coeffs = GetChebyDerCoeffs(coeffs);
ys_antider1_approx = Cheby2Phys(der1_coeffs(1:N+1));

% Take second antiderivative, enforcing similar boundary conditions
der2_coeffs = GetChebyDerCoeffs(der1_coeffs);
ys_antider2_approx = Cheby2Phys(der2_coeffs(1:N+1));

figure(1)
plot(xs,ys,'.-')
hold on
plot(xs,ys_approx,'o--')
hold off

figure(2)
plot(xs,ys_der1,'.-')
hold on
plot(xs,ys_antider1_approx,'o--')
hold off

figure(3)
plot(xs,ys_der2,'.-')
hold on
plot(xs,ys_antider2_approx,'o--')
hold off