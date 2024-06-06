function acd = GetChebyAntiDerCoeffs(coeffs,bound_value,xbound)
    % Finds antiderivative of given function in chebyshev space using its
    % coefficients
    % boundary condition used to evaluate constant of integration
    N = length(coeffs)-1;
    extCoeffs = [coeffs;0;0];
    ks = (1:N+1)';
    acd = [0;1./(2*ks).*(extCoeffs(1:end-2)-extCoeffs(3:end))];
    acd(2) = extCoeffs(1) - extCoeffs(3)/2;
    int_func = Cheby2Phys(acd(1:N+1));
    constant = bound_value - int_func(xbound);
    acd(1) = constant;
end