% modified from FT_CoRot_7.m

%%%%%%%%%%%%%%
% Parameters %
%%%%%%%%%%%%%%
DIR         = '/Volumes/DataSSD/OneDrive - UNSW/Hermite_Decomposition/ESO_HARPS/HD189733';
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
RV_FT = zeros(1,N_FILE);
h = figure; 
hold on
for i = 1:N_FILE
%     n = 1:size(FFT_frequency,2);
    n = (size(FFT_frequency,2)/2+1-10):(size(FFT_frequency,2)/2+1+10);
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

RV_FT = (RV_FT - mean(RV_FT)) * 1000;

% Time sequence %
h = figure; 
hold on
errorbar(MJD, RV_HARPS, RV_noise , 'r.', 'MarkerSize', 20)
errorbar(MJD, RV_FT, RV_noise , 'b.', 'MarkerSize', 20)
grid on
hold off
xlabel('Time [d]')
ylabel('RV [m/s]')
legend('HARPS', 'FT')
saveas(gcf,'5-Time_sequence','png')
% close(h)

% RV_HARPS vs RV_FT %
h = figure; 
plot(RV_HARPS, RV_FT, '.', 'MarkerSize', 10)
title('RV_{HARPS} vs RV_{FT}')
xlabel('RV_{HARPS} [m/s]')    
ylabel('RV_{FT} [m/s]')
saveas(gcf,'6-HARPS_vs_FT','png')
close(h)


% demo jitter % 
jitter_proto = (RV_HARPS - RV_FT') ;
t_smooth1 = linspace(54301.08, 54301.24, 1000);
y_smooth1 = FUNCTION_GAUSSIAN_SMOOTHING(MJD, jitter_proto, 1./RV_noise.^2, t_smooth1, 0.008);
t_smooth2 = linspace(54340.98, 54341.14, 1000);
y_smooth2 = FUNCTION_GAUSSIAN_SMOOTHING(MJD, jitter_proto, 1./RV_noise.^2, t_smooth2, 0.008);

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

    
    
% Fitting Part 1 %    
idx_t   = (MJD>54301.08) & (MJD<54301.24);
MJD1    = MJD(idx_t);    
weight1 = 1./RV_noise(idx_t).^2;
RV_HARPS1   = RV_HARPS(idx_t);
jitter1     = jitter_proto(idx_t);
jitter_smooth1  = FUNCTION_GAUSSIAN_SMOOTHING(MJD1, jitter1, weight1, MJD1, 0.010);

plot(MJD1, jitter1, 'o', MJD1, jitter_smooth1, '-')

MJD1        = MJD1 - min(MJD1);
idx_linear  = (MJD1>0.08);
MJD1_linear = MJD1(idx_linear);
RV_HARPS1_linear    = RV_HARPS1(idx_linear);
weight_linear       = weight1(idx_linear);
fun_b   = @(b) sum((b(1) + b(2) * MJD1_linear - RV_HARPS1_linear).^2 .*weight_linear);
b0  = [max(RV_HARPS1), -max(RV_HARPS1)/max(MJD1)];
b   = fminunc(fun_b,b0,options);

plot(MJD1, b(1) + b(2) * MJD1, '-', MJD1, RV_HARPS1, 'o')

% Jitter part %

RV_RM1  =  RV_HARPS1 - (b(1) + b(2) * MJD1);
fun_a   = @(a) sum((a(1) * jitter_smooth1 + a(2) - RV_RM1) .^2 .* weight1);
options = optimoptions(@fminunc,'Algorithm','quasi-newton');
a0 = [10, 0];
a = fminunc(fun_a,a0,options);

plot(MJD1, a(1) * jitter_smooth1 + a(2) + b(1) + b(2) * MJD1, '*', MJD1, RV_HARPS1, 'o')

res = a(1) * jitter_smooth1 + a(2) - RV_RM1;
plot(MJD1, res, '.')
    
%%%%%%%%%%%%%%%%%%
% Fitting Part 2 %   
%%%%%%%%%%%%%%%%%%

width   = 0.008;
idx_t   = (MJD>54340.98) & (MJD<54341.11);
MJD2    = MJD(idx_t);    
weight2 = RV_noise(idx_t);
RV_HARPS2   = RV_HARPS(idx_t);
RV_FT2      = RV_FT(idx_t);
RV_noise2   = RV_noise(idx_t);
jitter2     = jitter_proto(idx_t);

% Linear part %
MJD2    = MJD2 - min(MJD2);
idx_linear  = (MJD2 < 0.024) | (MJD2 > 0.107);
MJD2_linear = MJD2(idx_linear);
RV_HARPS2_linear    = RV_HARPS2(idx_linear);
weight_linear       = weight2(idx_linear);
fun_b = @(b) sum((b(1) + b(2) * MJD2_linear - RV_HARPS2_linear).^2 .* weight_linear);
b0 = [max(RV_HARPS2), -max(RV_HARPS2)/max(MJD2)];
b = fminunc(fun_b,b0,options);

h = figure;
hold on 
errorbar(MJD2, RV_HARPS2, RV_noise2, 'r.')
errorbar(MJD2, RV_FT2, RV_noise2, 'b*')
plot(MJD2, b(1) + b(2) * MJD2, '-')
xlim([-0.003 max(MJD2)+0.003])
legend('HARPS', 'FT')
title('RM effect as jitter')
xlabel('Time [d]')
ylabel('RV [m/s]')
saveas(gcf,'7-RM_effect1','png')
hold off
close(h)

% Jitter part %
t2_plot         = linspace(min(MJD2), max(MJD2), 1001);
jitter2_plot    = FUNCTION_GAUSSIAN_SMOOTHING(MJD2, jitter2, weight2, t2_plot, width);

h = figure;
hold on 
errorbar(MJD2, jitter2, RV_noise2, 'b*')
plot(t2_plot, jitter2_plot, 'b-')
xlim([-0.003 max(MJD2)+0.003])
xlabel('Time [d]')
ylabel('RV [m/s]')
legend('Proto jitter', 'Moving average')
saveas(gcf,'8-Proto_jitter2','png')
hold off
close(h)

jitter_smooth2  = FUNCTION_GAUSSIAN_SMOOTHING(MJD2, jitter2, weight2, MJD2, width);
RV_RM2 = RV_HARPS2 - (b(1) + b(2) * MJD2);
fun_a = @(a) sum((a(1) * jitter_smooth2 + a(2) - RV_RM2).^2 .* weight2);
options = optimoptions(@fminunc,'Algorithm','quasi-newton');
a0 = [10, 0];
a = fminunc(fun_a,a0,options);

h = figure;
subplot(2,1,1)       % add first plot in 2 x 1 grid
hold on 
errorbar(MJD2, RV_RM2, RV_noise2, 'r.')
errorbar(MJD2, a(1) * jitter_smooth2 + a(2), RV_noise2, 'b*')
ylabel('RV [m/s]')
title('RM effect as jitter')
legend('HARPS', 'FT')
hold off
subplot(2,1,2)       % add first plot in 2 x 1 grid
errorbar(MJD2, a(1) * jitter_smooth2 + a(2) - RV_RM2, RV_noise2, 'k.')
xlabel('Time [d]')
ylabel('Residual [m/s]')
saveas(gcf,'9-RM_fit2','png')
close(h)




    
dlmwrite('MJD_2012.txt', MJD-MJD(1))
dlmwrite('Jitter_model_2012.txt', jitter_proto)
dlmwrite('RV_HARPS_2012.txt', (RV_HARPS - mean(RV_HARPS)) * 1000)
dlmwrite('RV_FT_2012.txt', (RV_FT - mean(RV_FT))' * 1000)









