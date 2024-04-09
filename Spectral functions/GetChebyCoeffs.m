function coeffs = GetChebyCoeffs(ys)
    % For Chevyshev-Gauss-Lobatto domain, calculates coefficients for
    % Chebyshev expansion of data in ys
    N = length(ys)-1;
    c = [2,ones(1,N-1),2];
    coeffs = (2./(c*N).*sum((1./c') .* ys .* cos(pi*(0:N).*(0:N)'/(N))',1))';
end