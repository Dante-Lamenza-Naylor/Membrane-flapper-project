function smoothed_ys = GaussSmooth(ys,nondim_smoothScale)
    % Smooths membrane profile using a Gaussian kernel
    % - nondim_smoothScale:
    % Proportion of membrane chord length. Sets nodal "window size"
    % used in Gaussian kernel.

    N = length(ys);
    k = round(nondim_smoothScale*N);
    smoothed_ys =  smoothdata(ys,"gaussian",k);
    smoothed_ys([1,end]) = 0;
    
end