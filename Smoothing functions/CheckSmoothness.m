function [smooth,fft_xs,fft_ys] = CheckSmoothness(xs,ys,nondim_thres_contr, ...
                                      nondim_thres_freq)
% Checks smoothness of input vector ys on domain xs by taking fourier decomposition
% If high frequencies dominate: return false.
% - nondim_thres_contr:
% Proportion of maximum, in fourier decomposition of membrane
% profile's second spatial derivative in x, allowed for unwanted 
% frequencies.
% - nondim_thres_freq:
% Proportion of total sampling frequency, starting from 0 Hz, 
% permitted to have a relative contribution above 
% nondim_thres_contr.

% Y = nufft(X,t/T);
% n = length(t);
% f = (0:n-1)/n*Fs;
% plot(f,abs(Y))
% xlabel("f (Hz)")
% ylabel("abs(Y)(f)")

smooth = true;
N = length(xs);
xs = xs + xs(end); % To avoid negatives

L = xs(end) - xs(1); % Length of signal
Fs = N/L;            % Sampling frequency                                    
fft_ys = nufft(ys,xs*Fs);         % Compute fourier transform
% Output is reflected across midpoint; take first half for 
% relevant frequencies
fft_ys = abs(fft_ys(1:round(N/2)))/max(abs(fft_ys));
fft_xs = (0:round(N/2) - 1)'/N*Fs;
% Check for smoothness
threshold_frequency = nondim_thres_freq * Fs/2;
if max(fft_ys(fft_xs > threshold_frequency) > nondim_thres_contr)
    smooth = false;
end

    
end