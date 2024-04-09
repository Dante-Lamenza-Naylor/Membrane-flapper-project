function struct_force = RigidWingModel(y,y_targets,params)
    % Inputs should be column vectors
    
    k = params{8,2};
    struct_force = k*(y_targets-y);