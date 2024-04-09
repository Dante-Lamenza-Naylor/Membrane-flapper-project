function [dydx,dydx2,adydx] = ChebyCalculus(y)
    y = reshape(y,[],1);
    
    y_coeffs     = GetChebyCoeffs(y);
    dydx_coeffs  = GetChebyDerCoeffs(y_coeffs);
    dydx2_coeffs = GetChebyDerCoeffs(dydx_coeffs);
    adydx_coeffs = GetChebyAntiDerCoeffs(y_coeffs,0,1);
    adydx_coeffs = adydx_coeffs(1:end-1);
    
    dydx  = Cheby2Phys(dydx_coeffs);
    dydx2 = Cheby2Phys(dydx2_coeffs);
    adydx = Cheby2Phys(adydx_coeffs);