function cd = GetChebyDerCoeffs(coeffs)
    % Finds derivative of given function in chebyshev space using its
    % coefficients
    N = length(coeffs)-1;
    cd = zeros(N+2,1);
    cd(N) = N*coeffs(N);
    for k = (N+1:-1:3)
        cd(k-1) = cd(k+1) + 2*(k-1)*coeffs(k);
    end
    cd(1) = 0.5*cd(3) + coeffs(2);
    cd = cd(1:N+1);
end