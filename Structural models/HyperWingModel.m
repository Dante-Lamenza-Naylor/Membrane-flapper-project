function struct_force = HyperWingModel(y,y_targets,dydx,dydx2,params)
    % Inputs should be column vectors
    
    Ae       = params{2,2};
    lambda_0 = params{5,2};
    x        = params{10,2};

    kappa   =  dydx2.*(1 + (dydx).^2).^(-3/2);
    epsilon = lambda_0/2 * trapz(x,sqrt(1+(dydx).^2))-1;

    struct_force = 2*Ae*kappa*epsilon;