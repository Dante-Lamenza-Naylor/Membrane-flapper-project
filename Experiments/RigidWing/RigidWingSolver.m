function sol_return = RigidWingSolver(t_given,sol_vector,basket,params)

    N = params{1,2};
    R      = params{3,2};
    sigma  = params{4,2};
    amp    = params{6,2};
    k      = params{8,2};
    x      = params{10,2};

    basket.M_update_time(t_given)

    y = sol_vector(1:N+1);
    dydt = sol_vector(N+2:end);

    y_target_heave = amp * sin(sigma*t_given);
    y_targets = params{9,2} + y_target_heave;

    [dydx,dydx2,adydx] = ChebyCalculus(y);

    % struct_force = RigidWingModel(y,y_targets,params) - 1e2*4*R*k*dydt;
    struct_force = -amp*sigma^2*sin(sigma*t_given) * ones(N+1,1);
    dydt2 = struct_force; %./ (4*R);

    fluid_load = SpectralFluidModel1(x,dydt,dydt2,dydx,dydx2,params);

    sol_return = [dydt;dydt2];

    basket.STMsave(2,y);
    basket.STMsave(3,dydt);
    basket.STMsave(4,dydt2);
    basket.STMsave(5,struct_force);
    basket.STMsave(6,fluid_load);



    