% Use simulated spectra with planets AND stellar jitter: 
% /Volumes/DataSSD/SOAP_2/outputs/02.01/


%%%%%%%%%%%%%%
% Parameters %
%%%%%%%%%%%%%%
DIR         = '/Volumes/DataSSD/OneDrive - UNSW/Hermite_Decomposition/ESO_HARPS/LRa01_E2_0165';
file_list   = dir([DIR, '/4-ccf_dat/*.dat']);
file_name   = {file_list.name};
N_FILE      = size(file_name, 2);

MJD         = importdata([DIR, '/MJD.dat']);
RV_HARPS    = importdata([DIR, '/RV_HARPS.dat']);
x           = importdata([DIR, '/x.dat']);
weight      = importdata([DIR, '/weight.dat']);
grid_size   = 0.1;
Fs          = 1/grid_size;

if 0
% window function %
    window  = v0 * 0 + 1;
    bound   = 8;
    idx_w   = abs(v1) >= bound;
    window(idx_w)   = (cos((abs(v1(idx_w))-bound)/(10-bound)*pi) + 1) /2;

    h = figure;
    plot(v1, window)
    title('Window function')
    xlabel('Wavelength in RV [km/s]')
    ylabel('Window function')
    saveas(gcf,'0-Window_function','png')
    close(h)
    % window  = v1 * 0 + 1;
end




% estimate the size of array FFT_power
dat_name    = [DIR, '/4-ccf_dat/', char(file_name(1))];
A           = importdata(dat_name);
[aa, bb, yy]= FUNCTION_FFT(A, 0.1);
size1       = length(bb);
FFT_power   = zeros(size1, N_FILE);
Y           = zeros(size1, N_FILE);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stacked cross correlation function %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = figure;
hold on
for n = 1:N_FILE

    dat_name    = [DIR, '/4-ccf_dat/', char(file_name(n))];
    A           = importdata(dat_name);

%     A_spline    = A_spline .* window;
    plot(x, A)
    [FFT_frequency, FFT_power(:, n), Y(:, n)] = FUNCTION_FFT(A, Fs);
    
end     

hold off
title('Stacked cross correlation function')
xlabel('Wavelength in RV [km/s]')
ylabel('Normalized intensity')
saveas(gcf,'1-Line_Profile','png')
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
    xlabel('FT frequency (1 / velocity in wavelength)')
    ylabel('power')   
%     
    xlim([-0.5 0.5])
end 
hold off
saveas(gcf,'2-FT_power','png')
% saveas(gcf,'2-Differential_FT_power','png')
close(h)


%%%%%%%%%%%%%%%
% Phase angle %
%%%%%%%%%%%%%%%
h = figure;
plot(FFT_frequency, unwrap(angle(Y(:, 51))), '.')
title('Phase angle (Rotation phase = 0.51)')
xlabel('FT frequency (1 / velocity in wavelength)')
ylabel('Phase angle [radian]')
% xlim([-0.35 0.35])
saveas(gcf,'3-Phase_angle','png')
close(h)


% % phase 2D plot 
% [aa,bb,cc] = FUNCTION_FFT(A, Fs);
% h = figure;
% plot(real(cc), imag(cc), '.', real(Y(:,51)), imag(Y(:,51)), '.', 'MarkerSize', 10)
% legend('Noise free', 'SN=10000')
% grid on 
% xlim([-0.05 0.05])
% ylim([-0.005 0.005])
% xlabel('Real')    
% ylabel('Imaginary')
% title('Phase angle in complex plane')
% saveas(gcf,'7-Phase_angle_in_complex_plane_2','png')
% close(h)


%%%%%%%%%%%%%%%%%%%%%
% Phase angle -> RV %
%%%%%%%%%%%%%%%%%%%%%
slope = zeros(1,N_FILE);
RV_FT = zeros(1,N_FILE);
h = figure; 
hold on
for i = 1:N_FILE
%     n = 1:size(FFT_frequency,2);
    n = (size(FFT_frequency,2)/2+1-15):(size(FFT_frequency,2)/2+1+15);
    xx = FFT_frequency(n);
    yy = angle(Y(n, i)) - angle(Y(n, 1));
    plot(xx, yy, '-')
    title('Phase angle (relative to 1st epoch)')
    xlabel('FT frequency (1 / velocity in wavelength)')
    ylabel('Phase angle [radian]')
    
    % Phase angle -> RV
    p = polyfit(xx, unwrap(angle(Y(n, i)))',1);
    slope(i) = p(1);
    RV_FT1 = -slope(1) / (2*pi);
    RV_FT(i) = -slope(i) / (2*pi) - RV_FT1;
end
hold off
saveas(gcf,'4-Relative_phase_angle','png')
close(h)

if 1
    idx_obs2 = (MJD>5.55e4);
    MJD      = MJD(idx_obs2);
    RV_HARPS = RV_HARPS(idx_obs2);
    RV_FT    = RV_FT(idx_obs2);
    weight2   = weight(idx_obs2);
end

% MJD = MJD + 2400000.5 - 2450000;

% Time sequence %
h = figure; 
plot(MJD, (RV_HARPS - mean(RV_HARPS)) * 1000, 'b.', 'MarkerSize', 20); 
plot(MJD, (RV_FT - mean(RV_FT)) * 1000, 'r.', 'MarkerSize', 20)
rms_HARPS   = rms((RV_HARPS - mean(RV_HARPS)) * 1000);
rms_FT      = rms( (RV_FT - mean(RV_FT)) * 1000);
title('Time sequence')
xlabel('BJD')    
ylabel('RV [m/s]')
legend1 = ['rms_{HARPS} = ', num2str(rms_HARPS), ' m/s'];
legend2 = ['rms_{FT} = ', num2str(rms_FT), ' m/s'];
legend(legend1, legend2)
saveas(gcf,'5-Time_sequence','png')
close(h)

% RV_HARPS vs RV_FT %
h = figure; 
plot((RV_HARPS - mean(RV_HARPS)) * 1000, (RV_FT - mean(RV_FT)) * 1000, '.', 'MarkerSize', 10)
title('RV_{HARPS} vs RV_{FT}')
xlabel('RV_{HARPS} [m/s]')    
ylabel('RV_{FT} [m/s]')
saveas(gcf,'6-HARPS_vs_FT','png')
close(h)


% jitter % 
jitter_correction = ((RV_HARPS - mean(RV_HARPS)) * 1000 - (RV_FT - mean(RV_FT))' * 1000) ;
plot(MJD, jitter_correction, '.', 'MarkerSize', 20)
title('Jittere model')
xlabel('BJD')    
ylabel('Jittere model [m/s]')
dlmwrite('MJD_2012.txt', MJD-MJD(1))
dlmwrite('Jitter_model_2012.txt', jitter_correction)
dlmwrite('RV_HARPS_2012.txt', (RV_HARPS - mean(RV_HARPS)) * 1000)
dlmwrite('RV_FT_2012.txt', (RV_FT - mean(RV_FT))' * 1000)
dlmwrite('weight2.txt', weight2)





if 0
    plot(MJD, (RV_HARPS - mean(RV_HARPS)) * 1000, '.', MJD, jitter_correction, '.')
    plot(MJD, (RV_HARPS - mean(RV_HARPS)) * 1000, '*', MJD, (RV_HARPS - mean(RV_HARPS)) * 1000 - jitter_correction, 'o')
    legend('HARPS RV', 'RV+jitter')
    title('RV')
    xlabel('BJD')
    ylabel('RV [m/s]')

end



