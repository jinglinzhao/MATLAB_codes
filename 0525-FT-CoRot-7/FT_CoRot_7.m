% later modified from FT_HD189733.m


%%%%%%%%%%%%%%
% Parameters %
%%%%%%%%%%%%%%
DIR         = '/Volumes/DataSSD/OneDrive - UNSW/Hermite_Decomposition/ESO_HARPS/LRa01_E2_0165';
file_list   = dir([DIR, '/4-ccf_dat/*.dat']);
file_name   = {file_list.name};
N_FILE      = size(file_name, 2);

MJD         = importdata([DIR, '/MJD.dat']);
RV_HARPS    = importdata([DIR, '/RV_HARPS.dat']);                    % m/s
RV_HARPS    = (RV_HARPS - mean(RV_HARPS)) * 1000;
x           = importdata([DIR, '/x.dat']);
RV_noise      = importdata([DIR, '/RV_noise.dat']);
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


%%%%%%%%%%%%%%%%%%%%%
% Phase angle -> RV %
%%%%%%%%%%%%%%%%%%%%%
slope = zeros(1,N_FILE);
RV_FT = zeros(N_FILE, 1);
h = figure; 
    hold on
    for i = 1:N_FILE
        n = (size(FFT_frequency,2)/2+1-12):(size(FFT_frequency,2)/2+1+12);
        xx = FFT_frequency(n);
        yy = angle(Y(n, i)) - angle(Y(n, 1));
        plot(xx, yy, '-')
        title('Phase angle (relative to 1st epoch)')
        xlabel('FT frequency (1 / velocity in wavelength)')
        ylabel('Phase angle [radian]')

        % Phase angle -> RV
        p = polyfit(xx, unwrap(angle(Y(n, i)))',3);
        slope(i) = p(3);
        RV_FT1 = -slope(1) / (2*pi);
        RV_FT(i) = -slope(i) / (2*pi) - RV_FT1;
    end
    hold off
    saveas(gcf,'4-Relative_phase_angle','png')
close(h)

RV_FT = (RV_FT - mean(RV_FT)) * 1000;

% MJD = MJD + 2400000.5 - 2450000;

% Time sequence %
h = figure; 
    hold on
    errorbar(MJD, RV_HARPS, RV_noise , 'r.', 'MarkerSize', 20)
    errorbar(MJD, RV_FT, RV_noise , 'b.', 'MarkerSize', 20)
    grid on
    hold off
    title('Time sequence')
    xlabel('Time [d]')
    ylabel('RV [m/s]')
    legend1 = ['rms_{HARPS} = ', num2str(rms(RV_HARPS)), ' m/s'];
    legend2 = ['rms_{FT} = ', num2str(rms(RV_FT)), ' m/s'];
    legend(legend1, legend2)
    saveas(gcf,'5-Time_sequence','png')
close(h)

% RV_HARPS vs RV_FT %
h = figure; 
    plot(RV_HARPS, RV_FT, '.', 'MarkerSize', 10)
    title('RV_{HARPS} vs RV_{FT}')
    xlabel('RV_{HARPS} [m/s]')    
    ylabel('RV_{FT} [m/s]')
    saveas(gcf,'6-HARPS_vs_FT','png')
close(h)


% demo jitter % 
width   = 1.3;
jitter_proto = (RV_HARPS - RV_FT) ;

idx1    = (MJD<54900);
MJD1    = MJD(idx1);
t_smooth1 = linspace(min(MJD1), max(MJD1), 100001);
y_smooth1 = FUNCTION_GAUSSIAN_SMOOTHING(MJD, jitter_proto, 1./RV_noise.^2, t_smooth1, width);

idx2        = (MJD>54900);
MJD2        = MJD(idx2);
RV_HARPS2   = RV_HARPS(idx2);
RV_FT2      = RV_FT(idx2);
jitter2     = jitter_proto(idx2);
RV_noise2   = RV_noise(idx2);
t_smooth2   = linspace(min(MJD2), max(MJD2), 1001);
y_smooth2   = FUNCTION_GAUSSIAN_SMOOTHING(MJD, jitter_proto, 1./RV_noise.^2, t_smooth2, width);

h = figure; 
    hold on
    plot(t_smooth1, y_smooth1)
    plot(t_smooth2, y_smooth2)
    errorbar(MJD, jitter_proto, RV_noise, '.', 'MarkerSize', 20)
    hold off
    title('Jittere model')
    xlabel('BJD')    
    ylabel('Jittere model [m/s]')
    grid on
close(h)

dlmwrite('MJD_2012.txt', MJD2-min(MJD2))
dlmwrite('Jitter_model_2012.txt', jitter2)
dlmwrite('RV_HARPS_2012.txt', RV_HARPS2)
dlmwrite('RV_FT_2012.txt', RV_FT2)
dlmwrite('RV_noise_2012.txt', RV_noise2)


























