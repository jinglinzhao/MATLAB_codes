% Use Gaussian profile to study the line shift in FT

%%%%%%%%%%%%%%
% Parameters %
%%%%%%%%%%%%%%
SN          = 5000;
N_FILE      = 100;                               
grid_size    = 0.1;
Fs          = 1/grid_size;
v_max       = 10;
v_0         = (-v_max : grid_size : v_max)';          % km/s
% v              = (-100 : grid_size : 100)';      

% window function %
window  = v_0 * 0 + 1;
bound   = 8;
idx_w   = abs(v_0) >= bound;
window(idx_w)   = (cos((abs(v_0(idx_w))-bound)/(10-bound)*pi) + 1) /2;

h = figure;
plot(v_0, window)
title('Window function')
xlabel('Wavelength in RV [km/s]')
ylabel('Window function')
saveas(gcf,'0-Window_function','png')
close(h)

%%%%%%%%%%%%%%%%%%%
% Calculate Power %
%%%%%%%%%%%%%%%%%%%

% estimate the size of array FFT_power
A1          = exp(-v_0.^2/4^2);
[aa, bb, yy]= FUNCTION_FFT(A1, Fs);
size1       = length(bb);
FFT_power   = zeros(size1, N_FILE);
Y           = zeros(size1, N_FILE);
RV_noise    = zeros(1,N_FILE);
v_planet_array  = linspace(-3,3,101) / 1000.;
% v_planet_array  = 4 * sin(t/100.*1.8*2*pi + 1) * 0.001;
RV_gauss   = zeros(N_FILE,1);


h = figure;
hold on
for n = 1:N_FILE

    v_planet    = v_planet_array(n);
    v    = v_0 - v_planet;
    A    = exp(-v.^2/4^2);
    A    = A + normrnd(0, (1-A).^0.5/SN);  % add noise
%     A    = A .* window;
    plot(v, A)
    [FFT_frequency, FFT_power(:, n), Y(:, n)] = FUNCTION_FFT(A, Fs);
    
    f    = fit(v_0, A, 'a*exp(-((x-b)/c)^2)+d', 'StartPoint', [0.5, 0, 4, 0] );
    RV_gauss(n) = f.b * 1000;    
end     
hold off
title('Stacked Gaussian profile')
xlabel('x')
ylabel('y')
saveas(gcf,'1-Gaussian_Line_Profile','png')
close(h)


%%%%%%%%%%%%%%%%%%%%%%%%%%
% FT power in all epochs %
%%%%%%%%%%%%%%%%%%%%%%%%%%
h = figure; 
hold on
for n = 1:N_FILE
%     plot(FFT_frequency, FFT_power(:, n) - FFT_power(:, 51), '-')
%     title('Differential FT power in all epochs (overplot)')
%     xlabel('FT frequency (1 / velocity in wavelength)')
%     ylabel('Differential power') 
    
    plot(FFT_frequency, FFT_power(:, n), '.')
    title('FT power in all epochs (overplot)')
    xlabel('FT frequency')
    ylabel('power')   
%     
    xlim([-0.5 0.5])
end 
hold off
saveas(gcf,'2-FT_power','png')
close(h)


% Plot phase angle as a function of FT frequency (for zero RV shift)% 
h = figure;
hold on
% nn = 50:150;
plot(FFT_frequency, unwrap(angle(Y(:, 51))), 'r.')
% plot(FFT_frequency, angle(Y(:, 51)), 'b.')
hold off
title('Phase angle as a function of FT frequency (for zero RV shift)')
xlabel('FT frequency')
ylabel('Phase angle')    
% xlim([-0.35 0.35])
saveas(gcf,'3-Phase_angle','png')
close(h)


% Polyfit % 
slope = zeros(1,100);
h = figure; 
hold on
for i = 1:100
    n = (size(FFT_frequency,2)/2-0.5):(size(FFT_frequency,2)/2+1.5);
%     n = (size(FFT_frequency,2)/2-4.5):(size(FFT_frequency,2)/2+5.5);
%     n = 1:size(FFT_frequency,2);
    xx = FFT_frequency(n);ca
    yy = angle(Y(n, i)) - angle(Y(n, 51));
    yy(yy<-6)  = yy(yy<-6) + 2 * pi;
    yy(yy>6)  = yy(yy>6) - 2 * pi;
%     yy = unwrap(angle(Y(n, i))) - unwrap(angle(Y(n, 51)));
    plot(xx, yy, '-')
    title('Phase angle (relative to zero velocity)')
    xlabel('FT frequency')
    ylabel('Phase angle [rad]')   
    p = polyfit(xx,yy',1);
    slope(i) = p(1);
end
hold off
saveas(gcf,'4-Relative_phase_angle','png')
close(h)


% Compare with simulated RV % 
h = figure; 
xx = v_planet_array(1:100)*1000;
yy = -slope/(2*pi)*1000;
plot(xx, yy, 'o')
title('Recovered RV vs input RV')
xlabel('input RV (m/s)')    
ylabel('Recovered RV (m/s)')
saveas(gcf,'5-LINE_SHIFT_1','png')
close(h)

h = figure; 
xx          = v_planet_array(1:100)*1000;
yy          = -slope/(2*pi)*1000;
rms_gauss   = rms(yy - RV_gauss');
rms_FT      = rms(yy-xx);
plot(xx, yy - RV_gauss', '*', xx, yy-xx, 'o')
title('Recovered RV vs input RV')
xlabel('input RV (m/s)')    
ylabel('residual (m/s)')
legend1 = ['rms_{gauss} = ', num2str(rms_gauss), ' m/s'];
legend2 = ['rms_{FT} = ', num2str(rms_FT), ' m/s'];
legend(legend1, legend2)
saveas(gcf,'5-LINE_SHIFT_2','png')
close(h)


% More investigation in phase angle % 

h = figure; 
plot(FFT_frequency, angle(Y(:, 1)), '.', FFT_frequency, angle(Y(:, 51)), '.', FFT_frequency, angle(Y(:, 100)), '.')
legend('-3m/s', 'static', '2.9m/s shift')
xlabel('FT frequency')    
ylabel('Phase angle [rad]')
saveas(gcf,'3-Phase_angle_2','png')
close(h)


% phase 2D plot 
[aa,bb,cc] = FUNCTION_FFT(A1, Fs);
h = figure;
plot(real(cc), imag(cc), '.', real(Y(:,51)), imag(Y(:,51)), '.', 'MarkerSize', 10)
legend('Noise free', 'SN=10000')
grid on 
xlim([-0.05 0.05])
ylim([-0.005 0.005])
xlabel('Real')    
ylabel('Imaginary')
title('Phase angle in complex plane')
saveas(gcf,'7-Phase_angle_in_complex_plane_2','png')
close(h)





if 0
    h = figure; 
    Z = Y(:,51);
    Z(8:(201-6)) = 0;   % Set the high frequencies to zero
    ifft_Z = ifft(Z);   % inverse FT for only the low frequency components 
    plot(v_0, ifft_Z - A1, '.')     % plot the difference
    xlabel('x')
    ylabel('delta y')
    saveas(gcf,'6-Inverse_FT','png')
    close(h)

    figure; 
    hold on 
    for i =1:201
        plot(angle(Y(i, :)) - angle(Y(1, :)))
    end

end






