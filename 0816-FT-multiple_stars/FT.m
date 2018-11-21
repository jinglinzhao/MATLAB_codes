% modified from FT_HD189733.m

%%%%%%%%%%%%%%
% Parameters %
%%%%%%%%%%%%%%
star        = 'GJ699';
% star        = 'HD128621';
DIR         = ['/Volumes/DataSSD/OneDrive - UNSW/Hermite_Decomposition/ESO_HARPS/', star];
file_list   = dir([DIR, '/4-ccf_dat/*.dat']);
file_name   = {file_list.name};
N_FILE      = size(file_name, 2);

MJD         = importdata([DIR, '/MJD.dat']);
RV_HARPS    = importdata([DIR, '/RV_HARPS.dat']);                    % m/s
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
        plot(x, A-A1, 'k-')
    else
        plot(x, A-A1, 'r-')
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

% Determine the midpoint the equally divides the power spectrum %
cutoff_power= max(max(FFT_power)) * 0.001;
f_max       = max(FFT_frequency(FFT_power(:,1) > cutoff_power));
% f_max = 0.15;
n           = abs(FFT_frequency) <= f_max;
power_sum   = sum(FFT_power(n,1));
cum = 0;
for i = 0:fix(sum(n)/2)
    cum = cum + FFT_power(size(FFT_power,1)/2+i,1);
    if cum > power_sum/2
        break
    end
end
f_HL = FFT_frequency(size(FFT_power,1)/2+i);

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
close(h)

% Low-pass %
nl      = (FFT_frequency >= 0) & (FFT_frequency < f_HL);
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
% Time sequence %
h = figure; 
    hold on
%     errorbar(MJD, RV_HARPS, RV_noise , 'r.', 'MarkerSize', 20)
    errorbar(MJD, jitter_raw, RV_noise , 'b.', 'MarkerSize', 20)
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
    plot(RV_HARPS, RV_FTL, '.', 'MarkerSize', 10)
    title('RV_{HARPS} vs RV_{FT}')
    xlabel('RV_{HARPS} [m/s]')    
    ylabel('RV_{FT} [m/s]')
    saveas(gcf,'7-HARPS_vs_FT','png')
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

dlmwrite('jitter_smooth100.txt', y_smooth0)
jitter_raw = jitter_raw(idx);
dlmwrite('jitter_raw.txt', jitter_raw)

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
    % Fitting Part 1 %    
    idx_t   = (MJD>54301.08) & (MJD<54301.24);
    MJD1    = MJD(idx_t);    
    weight1 = 1./RV_noise(idx_t).^2;
    RV_HARPS1   = RV_HARPS(idx_t);
    RV_FT1      = RV_FT(idx_t);
    RV_noise1   = RV_noise(idx_t);
    jitter1     = jitter_proto(idx_t);

    % Linear part %
    MJD1        = MJD1 - min(MJD1);
    idx_linear  = (MJD1>0.08);
    MJD1_linear = MJD1(idx_linear);
    RV_HARPS1_linear    = RV_HARPS1(idx_linear);
    weight_linear       = weight1(idx_linear);
    fun_b   = @(b) sum((b(1) + b(2) * MJD1_linear - RV_HARPS1_linear).^2 .*weight_linear);
    options = optimoptions(@fminunc,'Algorithm','quasi-newton');
    b0  = [max(RV_HARPS1), -max(RV_HARPS1)/max(MJD1)];
    b   = fminunc(fun_b,b0,options);

    h = figure;
        hold on 
        errorbar(MJD1, RV_HARPS1, RV_noise1, 'r.')
        errorbar(MJD1, RV_FT1, RV_noise1, 'b*')
        plot(MJD1, b(1) + b(2) * MJD1, '-')
        xlim([-0.003 max(MJD1)+0.003])
        legend('HARPS', 'FT')
        title('Transit RV curve')
        xlabel('Time [d]')
        ylabel('RV [m/s]')
        saveas(gcf,'7-RM_effect1','png')
        hold off
    close(h)

    % Jitter part %

    t1_plot         = linspace(min(MJD1), max(MJD1), 1001);
    jitter1_plot    = FUNCTION_GAUSSIAN_SMOOTHING(MJD1, jitter1, weight1, t1_plot, width);

    h = figure;
        hold on 
        errorbar(MJD1, jitter1, RV_noise1, 'b*')
        plot(t1_plot, jitter1_plot, 'b-')
        xlim([-0.003 max(MJD1)+0.003])
        xlabel('Time [d]')
        ylabel('RV [m/s]')
        legend('Proto jitter', 'Moving average')
        title('Proto jitter')
        saveas(gcf,'8-Proto_jitter1','png')
        hold off
    close(h)

    jitter_smooth1  = FUNCTION_GAUSSIAN_SMOOTHING(MJD1, jitter1, weight1, MJD1, width);
    RV_RM1  =  RV_HARPS1 - (b(1) + b(2) * MJD1);
    fun_a   = @(a) sum((a(1) * jitter_smooth1 + a(2) - RV_RM1) .^2 .* weight1);
    options = optimoptions(@fminunc,'Algorithm','quasi-newton');
    a0 = [10, 0];
    a = fminunc(fun_a,a0,options);

    h = figure;
        subplot(2,1,1)       % add first plot in 2 x 1 grid
        hold on 
        errorbar(MJD1, RV_RM1, RV_noise1, 'r.')
        errorbar(MJD1, a(1) * jitter_smooth1 + a(2), RV_noise1, 'b*')
        xlim([-0.003 max(MJD1)+0.003])
        ylabel('RV [m/s]')
        title('Jitter recovery')
        legend('HARPS', 'FT')
        hold off
        subplot(2,1,2)       % add first plot in 2 x 1 grid
        errorbar(MJD1, a(1) * jitter_smooth1 + a(2) - RV_RM1, RV_noise1, 'k.')
        xlabel('Time [d]')
        ylabel('Residual [m/s]')
        xlim([-0.003 max(MJD1)+0.003])
        saveas(gcf,'9-RM_fit1','png')
    close(h)

    %%%%%%%%%%%%%%%%%%
    % Fitting Part 2 %   
    %%%%%%%%%%%%%%%%%%

    idx_t   = (MJD>54340.98) & (MJD<54341.11);
    MJD2    = MJD(idx_t);    
    weight2 = 1./RV_noise(idx_t).^2;
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
        subplot(3,1,1:2)
        hold on 
        scatter(MJD2, RV_FT2, 30, 'kD', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.5);
        scatter(MJD2, RV_HARPS2, 30, 'bo', 'MarkerFaceColor', 'b', 'MarkerFaceAlpha', 0.5);
    %     errorbar(MJD2, RV_FT2, RV_noise2, 'k.', 'MarkerSize', 0.1);
    %     errorbar(MJD2, RV_HARPS2, RV_noise2, 'b.', 'MarkerSize', 0.1);
        plot1 = plot(MJD2, b(1) + b(2) * MJD2, 'b--', 'LineWidth', 2);
        plot1.Color(4) = 0.4;        

        hold off
        xlim([-0.003 max(MJD2)+0.003])
        legend('FT', 'Gaussian')
    %     title('Transit RV curve')
        ylabel('RV [m/s]')
        set(gca,'fontsize', 15)
        set(gca,'xticklabel',[])    
        grid on

        % Jitter part %
        subplot(3,1,3)
        t2_plot         = linspace(min(MJD2), max(MJD2), 1001);
        jitter2_plot    = FUNCTION_GAUSSIAN_SMOOTHING(MJD2, jitter2, weight2, t2_plot, width);
        hold on
        scatter(MJD2, jitter2, 30, 'kD', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.5);
        errorbar(MJD2, jitter2, RV_noise2 * sqrt(2), 'k.', 'MarkerSize', 0.1)
        plot2 = plot(t2_plot, jitter2_plot, 'k-', 'LineWidth',2);
        plot2.Color(4) = 0.4;       
        hold off
        xlim([-0.003 max(MJD2)+0.003])
        xlabel('Time [d]')
        ylabel('\Delta RV [m/s]')
        set(gca,'fontsize', 15)
        saveas(gcf,'8-Proto_jitter2','png')
    close(h)

    jitter_smooth2  = FUNCTION_GAUSSIAN_SMOOTHING(MJD2, jitter2, weight2, MJD2, width);
    RV_RM2 = RV_HARPS2 - (b(1) + b(2) * MJD2);
    fun_a = @(a) sum((a(1) * jitter_smooth2 + a(2) - RV_RM2).^2 .* weight2);
    options = optimoptions(@fminunc,'Algorithm','quasi-newton');
    a0 = [10, 0];
    a = fminunc(fun_a,a0,options);

    h = figure;
    subplot(3,1,1:2)       % add first plot in 2 x 1 grid
        hold on 
        scatter(MJD2, a(1) * jitter2 + a(2), 30, 'ks', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.5);
        plot3 = plot(t2_plot, a(1) * jitter2_plot + a(2) , 'k', 'LineWidth', 2);
        plot3.Color(4) = 0.4;  
        scatter(MJD2, RV_RM2, 30, 'bo', 'MarkerFaceColor', 'b', 'MarkerFaceAlpha', 0.5);    
        errorbar(MJD2, a(1) * jitter2 + a(2), (RV_noise2 * sqrt(2) * a(1)), 'k.', 'MarkerSize', 0.1)
        errorbar(MJD2, RV_RM2, RV_noise2, 'b.', 'MarkerSize', 0.1)
        xlim([-0.003 max(MJD2)+0.003])
        ylabel('RV [m/s]')
        set(gca,'fontsize', 15)
        set(gca,'xticklabel',[])
        grid on 
        legend('Model', 'Smoothed model', 'RM effect (Jitter)')
        hold off

    subplot(3,1,3)       % add first plot in 2 x 1 grid
        hold on
    %     scatter(MJD2, a(1) * jitter2 + a(2) - RV_RM2, 30, 'ks', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.5)
    %     errorbar(MJD2, a(1) * jitter2 + a(2) - RV_RM2, (RV_noise2 * sqrt(2) * a(1)), 'k.', 'MarkerSize', 0.1)
        scatter(MJD2, a(1) * jitter_smooth2 + a(2) - RV_RM2, 30, 'bo', 'MarkerFaceColor', 'b', 'MarkerFaceAlpha', 0.5)
        plot5 = plot(MJD2, a(1) * jitter_smooth2 + a(2) - RV_RM2, 'k', 'LineWidth', 2);
        plot5.Color(4) = 0.4;
    %     plot4 = plot(MJD2, a(1) * jitter_smooth2 + a(2) - RV_RM2 , 'k', 'LineWidth', 2);
    %     plot4.Color(4) = 0.4;
        hold off
        xlabel('Time [d]')
        ylabel('Residual [m/s]')
    %     ylim([-30 30])
        xlim([-0.003 max(MJD2)+0.003])
        set(gca,'fontsize', 15)
        saveas(gcf,'9-RM_fit2','png')
    close(h)



    dlmwrite('MJD_2012.txt', MJD-MJD(1))
    dlmwrite('Jitter_model_2012.txt', jitter_proto)
    dlmwrite('RV_HARPS_2012.txt', (RV_HARPS - mean(RV_HARPS)) * 1000)
    dlmwrite('RV_FT_2012.txt', (RV_FT - mean(RV_FT))' * 1000)
end









