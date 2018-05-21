% Use Gaussian profile to study the line shift in FT

%%%%%%%%%%%%%%
% Parameters %
%%%%%%%%%%%%%%
SN              = 10000;
N_FILE          = 100;                               
grid_size       = 0.1;
Fs              = 1/grid_size;
% v              = (-10 : grid_size : 10)';          % km/s
v              = (-100 : grid_size : 100)';      

%%%%%%%%%%%%%%%%%%%
% Calculate Power %
%%%%%%%%%%%%%%%%%%%

% estimate the size of array FFT_power
A1          = exp(-v.^2/5^2);
[aa, bb, yy]= FUNCTION_FFT(A1, Fs);
size1       = length(bb);
FFT_power   = zeros(size1, N_FILE);
Y           = zeros(size1, N_FILE);
RV_noise    = zeros(1,N_FILE);
v_planet_array  = linspace(-3,3,101) / 1000.;


figure;
hold on
for n = 1:N_FILE

    % v_planet    = 3 * sin(n/100.*7*2*pi + 1) * 0.001;        % km/s
    v_planet    = v_planet_array(n);
    A           = exp(-v.^2/5^2);
%     plot(v0(idx),A(idx)-A1, '.')
    A_spline    = spline(v, A, v+v_planet);
%     A_spline    = A_spline - 0.0183;
    plot(v,A_spline)
    % A_spline    = A_spline + normrnd(0, (1-A_spline).^0.5/SN);  % add noise
    [FFT_frequency, FFT_power(:, n), Y(:, n)] = FUNCTION_FFT(A_spline, Fs);
end     
hold off


% FT power in all epochs % 
figure; 
hold on
for n = 1:N_FILE
    plot(FFT_frequency, FFT_power(:, n) - FFT_power(:, 51), '-')
    title('Differential FT power in all epochs (overplot)')
    xlabel('FT frequency (1 / velocity in wavelength)')
    ylabel('Differential power') 
    
%     plot(FFT_frequency, FFT_power(:, n))
%     title('FT power in all epochs (overplot)')
%     xlabel('FT frequency (1 / velocity in wavelength)')
%     ylabel('power')   
%     
    xlim([-0.2 0.2])
end 
hold off


% Plot phase angle as a function of FT frequency (for zero RV shift)% 
figure;
plot(FFT_frequency, angle(Y(:, 51)), '-')
title('Phase angle as a function of FT frequency (for zero RV shift)')
xlabel('FT frequency (1 / velocity in wavelength)')
ylabel('Phase angle')    
xlim([-0.35 0.35])


% Polyfit % 
slope = zeros(1,100);
figure; 
hold on
for i = 1:100
%     n = (size(FFT_frequency,2)/2-15):(size(FFT_frequency,2)/2+15);
%     n = (1025-40):(1025+40);    % plot for a particular frequency
    n = (size(FFT_frequency,2)/2-700):(size(FFT_frequency,2)/2+700);
    xx = FFT_frequency(n);
    yy = angle(Y(n, i)) - angle(Y(n, 51));
    plot(xx, yy, '-')
    title('Phase angle (relative to zero velocity)')
    xlabel('FT frequency (1 / velocity in wavelength)')
    ylabel('Phase angle')   
    p = polyfit(xx,yy',1);
    slope(i) = p(1);
end
hold off


% Compare with simulated RV % 
figure; 
xx = v_planet_array(1:100)*1000;
yy = slope/(2*pi)*1000;
plot(xx, yy, 'o')
title('Recovered RV vs input RV')
xlabel('input RV (m/s)')    
ylabel('Recovered RV (m/s)')

figure; 
xx = v_planet_array(1:100)*1000;
yy = slope/(2*pi)*1000;
plot(xx, yy-xx, 'o')
title('Recovered RV vs input RV')
xlabel('input RV (m/s)')    
ylabel('residual (m/s)')
