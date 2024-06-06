function sol_return = HeaveVacuumSolver(t_given,sol_vector,basket,params)

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

    % Boundary conditions
    y([1,end]) = y_targets([1,end]);

    % Do some derivatives with smoothing thrown in
    peak_threshold = 0.1;
    freq_threshold = 0.9;
    smooths = 0;
    while 1
        % [dydx,dydx2,adydx] = ChebyCalculus(y); UNSTABLE D:

        dydx = gradient(y)./gradient(x);
        dydx2 = gradient(dydx)./gradient(x);  % STABLE :D

        % Enforce edges
        dydx([1,end]) = 0;
        dydx2([1,end]) = 0;

        % Check smoothness
        [smooth,fft_xs,fft_ys] = CheckSmoothness(x,y,peak_threshold,freq_threshold); 
        if ~smooth
            smooths = smooths + 1;
            y = GaussSmooth(y,5);
            y([1,end]) = y_targets([1,end]);
        else
            break
        end
    end
    % disp(['smooths: ',num2str(smooths)])

    % ~~~~~~~~ SWAP WING MODELS HERE ~~~~~~~~~~~~~~~~~~~~~~
    struct_force = HyperWingModel(y,y_targets,dydx,dydx2,params);
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    dydt2 = struct_force ./ (4*R) - 0.1*dydt;

    sol_return = [dydt;dydt2];

    basket.STMsave(2,y);
    basket.STMsave(3,dydt);
    basket.STMsave(4,dydt2);
    basket.STMsave(5,struct_force);
    basket.STMsave(6,dydx2);
    basket.STMsave(7,fft_xs);
    basket.STMsave(8,fft_ys);



    