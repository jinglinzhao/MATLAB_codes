% modified from FT_HD189733.m

%%%%%%%%%%%%%%
% Parameters %
%%%%%%%%%%%%%%
% star        = 'GJ699';
% star        = 'HD224789';
star        = 'HD103720';
% star        = 'HD36051';
% star        = 'HD200143';
% star        = 'BD-213153';
% star        = 'HD216770';
% star        = 'HD189733';
% star        = 'HD22049';
% star        = 'HD128621';

DIR         = ['/Volumes/DataSSD/OneDrive - UNSW/Hermite_Decomposition/ESO_HARPS/', star];
file_list   = dir([DIR, '/4-ccf_dat/*.dat']);
file_name   = {file_list.name};
N_FILE      = size(file_name, 2);

MJD         = importdata([DIR, '/MJD.dat']);
RV_HARPS    = importdata([DIR, '/RV_HARPS.dat']);                  
RV_HARPS    = (RV_HARPS - mean(RV_HARPS)) * 1000;
x           = importdata([DIR, '/x.dat']);
RV_noise    = importdata([DIR, '/RV_noise.dat']);
grid_size   = 0.1;
Fs          = 1/grid_size;

% estimate the size of array FFT_power
dat_name    = [DIR, '/4-ccf_dat/', char(file_name(1))];
A1           = importdata(dat_name);
[aa, bb, yy]= FUNCTION_FFT(A1, 0.1);
size1       = length(bb);
FFT_power   = zeros(size1, N_FILE);
Y           = zeros(size1, N_FILE);


% construct a master template
h = figure;
hold on
A_tpl = zeros(size(A1));
for n = 1:N_FILE
    dat_name    = [DIR, '/4-ccf_dat/', char(file_name(n))];
    A           = importdata(dat_name);
    A_spline    = spline(x, A, x+(RV_HARPS(n)-RV_HARPS(1))/1000);
    A_tpl       = A_tpl + A_spline / RV_noise(n)^2;    
end     
A_tpl = A_tpl / sum(1./RV_noise.^2);
plot(x, A_tpl)
close(h)



h = figure;
hold on
array = 1:N_FILE;
for n = 1:N_FILE

    dat_name    = [DIR, '/4-ccf_dat/', char(file_name(n))];
    A           = importdata(dat_name);
%     A_spline    = spline(x, A, x+(BI(n)-BI(1))/1000);
    A_spline    = spline(x, A, x+(RV_HARPS(n)-RV_HARPS(1))/1000);
    plot(x, A_spline - A_tpl, 'k-')
end     
errorbar(median(x), 0, 1/7676, 'r', 'LineWidth',3')

hold off
% title('Stacked cross correlation function')
% ylim([-0.001 0.1])
set(gca,'fontsize',24)
xlabel('km/s')
ylabel('Normalized intensity')
% title('Line profile (stacked)')
title('Centred differential line profile')
close(h)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stacked cross correlation function %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = figure;
hold on
for n = 1:N_FILE

    dat_name    = [DIR, '/4-ccf_dat/', char(file_name(n))];
    A           = importdata(dat_name);

%     A_spline    = A_spline .* window;
%     if mod(n,5) == 1
%         plot(x, A-A1, 'k-')
%     end

    if MJD(n)<57161
        plot(x, A-A_tpl, 'k-')
    else
        plot(x, A-A_tpl, 'r-')
    end
    
    if strcmp(star, 'HD103720')
        A = spline(x, A, x+(RV_HARPS(n)-RV_HARPS(1))/1000);
    end
    
    [FFT_frequency, FFT_power(:, n), Y(:, n)] = FUNCTION_FFT(A, Fs);
end     
hold off
% title('Stacked cross correlation function')
% ylim([-0.001 0.1])
set(gca,'fontsize',20)
xlabel('km/s')
ylabel('Normalized intensity')
title(star)
% saveas(gcf,'1-Line_Profile','png')
saveas(gcf,'1-Differential_line_Profile','png')
close(h)


%%%%%%%%%%%%
% Archived %
%%%%%%%%%%%%
% cutoff_power= max(max(FFT_power)) * 0.001;
% f_max       = max(FFT_frequency(FFT_power(:,1) > cutoff_power));
% n           = abs(FFT_frequency) <= f_max;
% power_sum   = sum(FFT_power(n,1));
% cum = 0;
% for i = 0:fix(sum(n)/2)
%     cum = cum + FFT_power(size(FFT_power,1)/2+i,1);
%     if cum > power_sum/2
%         break
%     end
% end
% f_HL = FFT_frequency(size(FFT_power,1)/2+i);


% Determine the midpoint the equally divides the power spectrum %
% cutoff_power= max(max(FFT_power)) * 0.005; % HD224789
cutoff_power= max(max(FFT_power)) * 0.001; % HD103720
f_max       = max(FFT_frequency(FFT_power(:,1) > cutoff_power));
n           = abs(FFT_frequency) <= f_max;
power_sum   = sum(FFT_power(n,1));
% half power %
cum = 0;
for i = 1:fix(sum(n)/2)
    cum = cum + FFT_power(size(FFT_power,1)/2+1+i,1);
    if cum > power_sum/4
        break
    end
end
f_HL = FFT_frequency(size(FFT_power,1)/2+1+i);
% f_HL = f_max * 0.5;
% n1 = abs(FFT_frequency) <= f_HL;
% power_sum1   = sum(FFT_power(n1,1)); % 95% of power_sum



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
    
    plot(FFT_frequency, FFT_power(:, n), 'k')
%     title('Stacked power spectrum')
end 
hold off
xlabel('\xi [s/km]')
ylabel('Power')   
xlim([-0.26 0.26])
set(gca,'fontsize',20)
saveas(gcf,'2-FT_power','png')
% saveas(gcf,'2-Differential_FT_power','png')
close(h)


%%%%%%%%%%%%%%%
% Phase angle %
%%%%%%%%%%%%%%%
h = figure;
plot(FFT_frequency, unwrap(angle(Y(:, 1))), '.')
title('Phase angle (Rotation phase = 0.51)')
xlabel('FT frequency (1 / velocity in wavelength)')
ylabel('Phase angle [radian]')
xlim([-0.35 0.35])
saveas(gcf,'3-Phase_angle','png')
close(h)


%%%%%%%%%%%%%%%%%%%%%
% Phase angle -> RV %
%%%%%%%%%%%%%%%%%%%%%
n       = abs(FFT_frequency) <= f_max;
slope = zeros(1,N_FILE);
RV_FT  = zeros(1,N_FILE);
wegihted_velocity = zeros(1,N_FILE);
h = figure; 
hold on
for i = 1:N_FILE
    xx  = FFT_frequency(n);
    yy  = angle(Y(n, i)) - angle(Y(n, 1));
%     if mod(i,5) == 1
        plot(xx, yy, 'k-')
%     end
    % Phase angle -> RV
    weight = FFT_power(n,i)';
    [fitresult, gof] = createFit(xx, yy', weight);
    slope(i) = fitresult.p1;
    RV_FT(i) = -slope(i) / (2*pi);
end
hold off
set(gca,'fontsize',20)
xlabel('\xi [s/km]')
ylabel('\Delta \phi [radian]')
saveas(gcf,'4-Relative_phase_angle','png')


% Low-pass %
nl      = (FFT_frequency >= 0) & (FFT_frequency <= f_HL);
RV_FTL  = zeros(1,N_FILE);
h       = figure; 
hold on
for i = 1:N_FILE
    xx  = FFT_frequency(nl);
    yy  = angle(Y(nl, i)) - angle(Y(nl, 1));
%     if mod(i,5) == 1
        plot(xx, yy, 'k-')
%     end
    % Phase angle -> RV
    weight = FFT_power(nl,i)';
    [fitresult, gof] = createFit(xx, yy', weight);
    slope(i) = fitresult.p1;
    RV_FTL(i) = -slope(i) / (2*pi);
end
hold off
set(gca,'fontsize',20)
xlabel('\xi [s/km]')
ylabel('\Delta \phi [radian]')
title('Low-pass')
saveas(gcf,'4-Relative_phase_angle_L','png')
close(h)

% high-pass % 
n       = (FFT_frequency >= f_HL) & (FFT_frequency <= f_max);
RV_FTH  = zeros(1,N_FILE);
h       = figure; 
hold on
for i = 1:N_FILE
    xx  = FFT_frequency(n);
    yy  = angle(Y(n, i)) - angle(Y(n, 1));
%     if mod(i,5) == 1
        plot(xx, yy, 'k-')
%     end
    % Phase angle -> RV
    weight = FFT_power(n,i)';
    [fitresult, gof] = createFit(xx, yy', weight);
    slope(i) = fitresult.p1;
    RV_FTH(i) = -slope(i) / (2*pi);    
end
hold off
set(gca,'fontsize',20)
xlabel('\xi [s/km]')
ylabel('\Delta \phi [radian]')
title('High-pass')
saveas(gcf,'4-Relative_phase_angle_H','png')
close(h)

% RV_FT   = RV_FT' * 1000;
% RV_FTL  = RV_FTL' * 1000;
% RV_FTH  = RV_FTH' * 1000;

RV_FT   = (RV_FT - mean(RV_FT))' * 1000;
RV_FTL  = (RV_FTL - mean(RV_FTL))' * 1000;
RV_FTH  = (RV_FTH - mean(RV_FTH))' * 1000;

dlmwrite('GG.txt', RV_HARPS)
dlmwrite('XX.txt', RV_FT)
dlmwrite('YY.txt', RV_FTL)
dlmwrite('ZZ.txt', RV_FTH)

jitter_raw  = RV_HARPS-RV_FTL;
% jitter_raw  = RV_FTH - RV_HARPS;
% MJD = MJD - 50000;

% visual check
h = figure; 
    hold on
    errorbar(MJD, RV_HARPS, RV_noise , 'k.', 'MarkerSize', 20)
    errorbar(MJD, RV_FT, RV_noise , 'r.', 'MarkerSize', 20)
    grid on
    hold off
    xlabel('Time [d]')
    ylabel('RV [m/s]')
    title(star)
close(h)

% Time sequence %
h = figure; 
    hold on
    errorbar(MJD, RV_HARPS, RV_noise , 'k.', 'MarkerSize', 20)
    errorbar(MJD, jitter_raw, RV_noise , 'b.', 'MarkerSize', 20)
%     errorbar(MJD, RV_FTH, RV_noise , 'r.', 'MarkerSize', 20)
    grid on
    hold off
    xlabel('Time [d]')
    ylabel('RV [m/s]')
%     legend('HARPS', 'FT')
    title(star)
    saveas(gcf,'5-Time_sequence','png')
close(h)

% proto-jitter %
h = figure; 
    plot(RV_HARPS-RV_FTL, RV_FTH-RV_HARPS, '.', 'MarkerSize', 10)
    xlabel('RV_{HARPS} - RV_{FT,L} [m/s]')    
    ylabel('RV_{FT,H} - RV_{HARPS} [m/s]')
    saveas(gcf,'6-proto-jitter','png')
close(h)

% RV_HARPS vs RV_FT %
h = figure; 
    hold on
    plot(RV_HARPS, RV_FTL, 'b.', 'MarkerSize', 10)
%     plot(RV_HARPS, RV_FTH, 'r.', 'MarkerSize', 10)
    hold off
    title('RV_{HARPS} vs RV_{FT}')
    xlabel('RV_{HARPS} [m/s]')    
    ylabel('RV_{FT} [m/s]')
    saveas(gcf,'7-HARPS_vs_FT','png')
close(h)

% RV_HARPS vs RV_FT %
h = figure; 
    hold on
    plot(RV_HARPS, jitter_raw, 'b.', 'MarkerSize', 10)
    hold off
    title('RV_{HARPS} vs \DeltaRV_{FT,L}')
    xlabel('RV_{HARPS} [m/s]')    
    ylabel('\DeltaRV_{FT,L} [m/s]')
close(h)


%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(star,'HD85390')
%%%%%%%%%%%%%%%%%%%%%%%%%    
    width   = 1;
    idx     = MJD < 57300;
    jitter_proto = RV_FTH(idx) - RV_FTL(idx);
    MJD1    = MJD(idx);
    y_smooth0 = FUNCTION_GAUSSIAN_SMOOTHING(MJD1, jitter_proto, 1./RV_noise(idx).^2, MJD1, width);
    t_smooth1 = linspace(min(MJD1), max(MJD1), 100000);
    y_smooth1 = FUNCTION_GAUSSIAN_SMOOTHING(MJD1, jitter_proto, 1./RV_noise(idx).^2, t_smooth1, width);

    h = figure; 
        hold on
        p7 = plot(t_smooth1, y_smooth1, 'LineWidth', 2); p7.Color(4)=0.4;
        scatter(MJD1, y_smooth0, 15, 'ro', 'MarkerFaceColor', 'r', 'MarkerFaceAlpha', 0.1)
        errorbar(MJD1, jitter_proto, 2*RV_noise(idx), 'b.', 'MarkerSize', 20)
        hold off
        title('Jittere model')
        xlabel('shifted BJD')    
        ylabel('Radial velocity [m/s]')
        grid on
        xlim([min(MJD1)-100 max(MJD1)+100])
    %     pbaspect([3 1 1])
        set(gca,'fontsize',20)
        saveas(gcf,'7-Jitter_Model','png')
    close(h)
end

% dlmwrite('jitter_smooth100.txt', y_smooth0)
% jitter_raw = jitter_raw(idx);
% dlmwrite('jitter_raw.txt', jitter_raw)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(star,'LRa01_E2_0165')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    width   = 1.3;
    jitter_proto = jitter_raw;

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
end

if 0
    dlmwrite('MJD_2012.txt', MJD-MJD(1))
    dlmwrite('Jitter_model_2012.txt', jitter_proto)
    dlmwrite('RV_HARPS_2012.txt', (RV_HARPS - mean(RV_HARPS)) * 1000)
    dlmwrite('RV_FT_2012.txt', (RV_FT - mean(RV_FT))' * 1000)
end









