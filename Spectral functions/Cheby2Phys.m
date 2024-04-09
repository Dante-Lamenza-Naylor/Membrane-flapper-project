function ys_approx = Cheby2Phys(coeffs)
% Transforms function from chebyshev to physical space
% Domain must be Chebyshev-Gauss-Lobatto collocation points
N = length(coeffs)-1;
ys_approx = sum(coeffs.*cos(pi*(0:N).*(0:N)'/(N))')';