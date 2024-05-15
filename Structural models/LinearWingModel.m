function struct_force = LinearWingModel(y,y_targets,dydx,dydx2,params)
    % Inputs should be column vectors
    
    Ae = params{2,2};
    R  = params{3,2};
    lambda_0 = params{5,2};

    struct_force = 2*Ae*(lambda_0 - 1)*dydx2;