trial_func = @(x)x;
trial_func_der1 = @(x)ones(length(x),1);
trial_func_der2 = @(x)zeros(length(x),1);
trial_func_der3 = @(x)zeros(length(x),1);
trial_func_der4 = @(x)zeros(length(x),1);
N=128;
xs = cos(pi*(0:N)/(N)); % Must be Chebyshev-Gauss-Lobatto collocation points
ys = trial_func(xs)';
ys_der1 = trial_func_der1(xs)';
ys_der2 = trial_func_der2(xs)';
ys_der3 = trial_func_der3(xs)';
ys_der4 = trial_func_der4(xs)';

% Convert to spectral space
coeffs = GetChebyCoeffs(ys);
ys_approx = Cheby2Phys(coeffs);
ys_sin_approx = 2*sum(coeffs(2:end).*sin(pi*(1:N).*(1:N)'/N)');
% ys_spline = makima(xs,ys);

% Take antiderivative, enforcing boundary condition =0 at x=-1
der1_test_coeffs = GetChebyCoeffs(trial_func_der1(xs)');
der1_coeffs = GetChebyDerCoeffs(coeffs);
ys_antider1_approx = Cheby2Phys(der1_coeffs(1:N+1));
ys_antider1_findiff = gradient(ys)./gradient(xs');

% Take second antiderivative, enforcing similar boundary conditions
der2_coeffs = GetChebyDerCoeffs(der1_coeffs);
ys_antider2_approx = Cheby2Phys(der2_coeffs(1:N+1));
ys_antider2_findiff = gradient(ys_antider1_findiff)./gradient(xs');

der3_coeffs = GetChebyDerCoeffs(der2_coeffs);
ys_antider3_approx = Cheby2Phys(der3_coeffs(1:N+1));
ys_antider3_findiff = gradient(ys_antider2_findiff)./gradient(xs');

der4_coeffs = GetChebyDerCoeffs(der3_coeffs);
ys_antider4_approx = Cheby2Phys(der4_coeffs(1:N+1));
ys_antider4_findiff = gradient(ys_antider3_findiff)./gradient(xs');

figure(1)
plot(xs,ys,'.-')
hold on
plot(xs,ys_approx,'o--')
hold off
title("0th derivative")

figure(2)
plot(xs,ys_der1,'.-')
hold on
plot(xs,ys_antider1_approx,'o--')
plot(xs,ys_antider1_findiff,'x--')
hold off
title("1st derivative")
legend('ana','cheb','findiff')

figure(3)
plot(xs,ys_der2,'.-')
hold on
plot(xs,ys_antider2_approx,'o--')
plot(xs,ys_antider2_findiff,'x--')
hold off
title("2nd derivative")
legend('ana','cheb','findiff')

figure(4)
plot(xs,ys_der3,'.-')
hold on
plot(xs,ys_antider3_approx,'o--')
plot(xs,ys_antider3_findiff,'x--')
hold off
title("3rd derivative")
legend('ana','cheb','findiff')

figure(5)
plot(xs,ys_der4,'.-')
hold on
plot(xs,ys_antider4_approx,'o--')
plot(xs,ys_antider4_findiff,'x--')
hold off
title("4th derivative")
legend('ana','cheb','findiff')

figure(6)
semilogy(xs,abs(ys-ys_approx)./ys,'o-')
title("0th der error")

figure(7)
semilogy(xs,abs(ys_der1-ys_antider1_approx)./abs(ys_der1),'o--')
hold on
semilogy(xs,abs(ys_der1-ys_antider1_findiff)./abs(ys_der1),'x--')
hold off
title("1st der error")
legend('cheb','findiff')

figure(8)
semilogy(xs,abs(ys_der2-ys_antider2_approx)./abs(ys_der2),'o--')
hold on
semilogy(xs,abs(ys_der2-ys_antider2_findiff)./abs(ys_der2),'x--')
hold off
title("2nd der error")
legend('cheb','findiff')

figure(9)
semilogy(xs,abs(ys_der3-ys_antider3_approx)./abs(ys_der3),'o--')
hold on
semilogy(xs,abs(ys_der3-ys_antider3_findiff)./abs(ys_der3),'x--')
hold off
title("3rd der error")
legend('cheb','findiff')

figure(10)
semilogy(xs,abs(ys_der4-ys_antider4_approx)./abs(ys_der4),'o--')
hold on
semilogy(xs,abs(ys_der4-ys_antider4_findiff)./abs(ys_der4),'x--')
hold off
title("4th der error")
legend('cheb','findiff')
