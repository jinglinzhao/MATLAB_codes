function [fshift,power2] = FUNCTION_FFT_iFFT(X, Fs)

% To use the fft function to convert the signal to the frequency domain, 
% first identify a new input length that is the next power of 2 from the 
% original signal length. This will pad the signal X with trailing zeros 
% in order to improve the performance of fft.

n       = 2^nextpow2(length(X))*8;
fshift  = Fs*((0:n-1)/n-0.5);
% fshift  = transpose(fshift);

% Power is the squared magnitude of a signal's Fourier transform,
% normalized by the number of frequency samples.
Y = fft(X, n);
Yshift = fftshift(Y);    
power = abs(Yshift)/n; % squared? 

y       = ifft(power);
yshift  = fftshift(y);   
power2  = abs(yshift)*n; 
idx = (fshift<10) & (fshift>-10);
figure; plot(fshift(idx),power2(idx))

end

% synchronize the file to the data folder 
% rsync /Volumes/DataSSD/MATLAB_codes/Project180129-FT/FUNCTION_FFT_noise.m /Volumes/DataSSD/OneDrive\ -\ UNSW/Hermite_Decomposition/ESO_HARPS/code
% rsync /Volumes/DataSSD/MATLAB_codes/Project180129-FT/FUNCTION_FFT_noise.m /Volumes/DataSSD/MATLAB_codes/Project180131-FT_SOAP
% 
% 