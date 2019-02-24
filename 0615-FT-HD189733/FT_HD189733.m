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
A1           = importdata(dat_name);
[aa, bb, yy]= FUNCTION_FFT(A1, 0.1);
size1       = length(bb);
FFT_power   = zeros(size1, N_FILE);
Y           = zeros(size1, N_FILE);


%%%%%%%%
% Test %
%%%%%%%%
if 0
    n0 = (MJD2<0.025) | (MJD2>0.104);
    n1 = (MJD2<0.062) & (MJD2>0.025);
    n2 = (MJD2>0.062) & (MJD2<0.104);
    figure;
    hold on 
    plot(MJD2(n0), RV_HARPS2(n0)-BI(n0), 'k.')
    plot(MJD2(n1), RV_HARPS2(n1)-BI(n1), 'ro')
    plot(MJD2(n2), RV_HARPS2(n2)-BI(n2), 'b+')
    hold off 
end 

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


h = figure;
hold on
array = 1:N_FILE;
for n = 1:N_FILE

    dat_name    = [DIR, '/4-ccf_dat/', char(file_name(n))];
    A           = importdata(dat_name);
%     A_spline    = spline(x, A, x+(BI(n)-BI(1))/1000);
    A_spline    = spline(x, A, x+(RV_HARPS(n)-RV_HARPS(1))/1000);
    
%     if ismember(n, array(n0))
%         plot(x, A_spline - A_tpl, 'k-')
% %         
%     end
%     if ismember(n, array(n1))
%         plot(x, A_spline - A_tpl, 'b-')
% %         plot(x, A, 'b-')
%     end    
%     if ismember(n, array(n2))
%         plot(x, A_spline - A_tpl, 'r-')
% %         plot(x, A, 'r-')
%     end
    
    % NEW %
%     A_spline    = spline(v1, A_spline, v1+f.b);    
    [FFT_frequency, FFT_power(:, n), Y(:, n)] = FUNCTION_FFT(A_spline, Fs);
    
end     
hold off
set(gca,'fontsize', 20)
xlabel('km/s')
ylabel('Normalized intensity')
% title('Line profile of HD 189733')
% title({'Residual line profile', 'of HD 189733'})
title({'Centred residual line profile', 'of HD 189733'})
ylim([-0.0008 0.0008])
% ylim([0 0.09])
% saveas(gcf,'HD189733LineProfile','png')
saveas(gcf,'HD189733DLP2','png')



saveas(gcf,'Line_Profile','png')
% saveas(gcf,'1-Line_Profile','png')
saveas(gcf,'1-Differential_line_Profile','png')
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
        plot(x, A-A1, 'k-')
%     end
%     plot(x, A-A1, 'k')
    [FFT_frequency, FFT_power(:, n), Y(:, n)] = FUNCTION_FFT(A, Fs);
    
end     

hold off
% title('Stacked cross correlation function')
% ylim([-0.001 0.1])
set(gca,'fontsize',20)
xlabel('km/s')
ylabel('Normalized intensity')
% saveas(gcf,'1-Line_Profile','png')
saveas(gcf,'1-Differential_line_Profile','png')
close(h)

% Determine the midpoint the equally divides the power spectrum %
cutoff_power= max(max(FFT_power)) * 0.001;
f_max       = max(FFT_frequency(FFT_power(:,1) > cutoff_power));
n           = abs(FFT_frequency) <= f_max;
power_sum   = sum(FFT_power(n,1));
cum = 0;
for i = 1:fix(sum(n)/2)
    cum = cum + FFT_power(size(FFT_power,1)/2+1+i,1);
    if cum > power_sum/4
        break
    end
end
f_HL = FFT_frequency(size(FFT_power,1)/2+1+i);

if 0 % old 
    for i = 0:fix(sum(n)/2)
        cum = cum + FFT_power(size(FFT_power,1)/2+i,1);
        if cum > power_sum/2
            break
        end
    end
    f_HL = FFT_frequency(size(FFT_power,1)/2+i);
end

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
xlim([-0.15 0.1501])
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
% xlim([-0.35 0.35])
saveas(gcf,'3-Phase_angle','png')
close(h)


%%%%%%%%%%%%%%%%%%%%%
% Phase angle -> RV %
%%%%%%%%%%%%%%%%%%%%%
n       = abs(FFT_frequency) <= f_max;
slope   = zeros(1,N_FILE);
RV_FT   = zeros(1,N_FILE);
RV_FT_err  = zeros(1,N_FILE);
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
    ci              = confint(fitresult,0.95);
    RV_FT_err(i)    = abs(diff(ci(:,1))*1000 / (4*pi));    
end
hold off
set(gca,'fontsize',20)
xlabel('\xi [s/km]')
ylabel('\Delta \phi [radian]')
saveas(gcf,'4-Relative_phase_angle','png')
close(h)

% Low-pass %
nl      = (FFT_frequency >= 0) & (FFT_frequency <= f_HL);
RV_FTL  = zeros(1,N_FILE);
RV_FTL_err  = zeros(1,N_FILE);
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
    ci              = confint(fitresult,0.95);
    RV_FTL_err(i)   = abs(diff(ci(:,1))*1000 / (4*pi));        
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
RV_FTH_err  = zeros(1,N_FILE);
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
    ci              = confint(fitresult,0.95);
    RV_FTH_err(i)   = abs(diff(ci(:,1))*1000 / (4*pi));            
end
hold off
set(gca,'fontsize',20)
xlabel('\xi [s/km]')
ylabel('\Delta \phi [radian]')
title('High-pass')
saveas(gcf,'4-Relative_phase_angle_H','png')
close(h)

XX = (RV_FT-mean(RV_FT))'*1000;
YY = (RV_FTL-mean(RV_FTL))'*1000;
ZZ = (RV_FTH-mean(RV_FTH))'*1000;

idx_t   = (MJD>54340.98) & (MJD<54341.11);
% idx_t = (MJD > 53986) | (MJD < 53990);
GGG = RV_HARPS(idx_t);
XXX = XX(idx_t);
YYY = YY(idx_t);
ZZZ = ZZ(idx_t);

BI = b(1) + b(2) * MJD2;

dlmwrite('GG.txt', GGG)
dlmwrite('XX.txt', XXX)
dlmwrite('YY.txt', YYY)
dlmwrite('ZZ.txt', ZZZ)
dlmwrite('BI.txt', BI)
dlmwrite('RV_noise.txt', RV_noise(idx_t))

dlmwrite('XY.txt', jitter_smooth2)


 



% Time sequence %
h = figure; 
    hold on
%     errorbar(MJD-min(MJD), RV_HARPS, RV_noise , 'r.', 'MarkerSize', 20)
%     errorbar(MJD-min(MJD), YY, RV_noise , 'b.', 'MarkerSize', 20)
    errorbar(MJD-min(MJD), ZZ, RV_noise , 'b.', 'MarkerSize', 20)
    grid on
    hold off
    xlabel('Time [d]')
    ylabel('RV [m/s]')
    legend('HARPS', 'FT')
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
width   = 0.008;
% RV_FT = ZZ;  
% RV_FT = YY; 
jitter_proto = ZZZ;
% jitter_proto = (ZZ - RV_HARPS) + (RV_HARPS - YY); 
% jitter_proto = (ZZ - RV_HARPS);
% jitter_proto = (RV_HARPS - YY);

t_smooth1 = linspace(54301.08, 54301.24, 1000);
y_smooth1 = FUNCTION_GAUSSIAN_SMOOTHING(MJD, jitter_proto, 1./RV_noise.^2, t_smooth1, width);
% Section 4
t_smooth2 = linspace(54340.98, 54341.14, 1000);
y_smooth2 = FUNCTION_GAUSSIAN_SMOOTHING(MJD, jitter_proto, 1./RV_noise.^2, t_smooth2, width);
% Section 2
t_smooth  = linspace(53986, 53986.18, 1000);
y_smooth = FUNCTION_GAUSSIAN_SMOOTHING(MJD, jitter_proto, 1./RV_noise.^2, t_smooth, width);


h = figure; 
    hold on
%     plot(t_smooth1, y_smooth1)
%     plot(t_smooth2, y_smooth2)
%     plot(t_smooth, y_smooth)
    errorbar(MJD - min(MJD), jitter_proto, RV_noise*3^0.5, '.', 'MarkerSize', 20)
    hold off
    title('Jittere model')
    xlabel('BJD')    
    ylabel('Jittere model [m/s]')
    grid on
close(h)




%%%%%%%%%%%%%%%%%%
% Fitting Part 0 %   
%%%%%%%%%%%%%%%%%%
weight = 1./RV_noise.^2;
RV_HARPS   = RV_HARPS - mean(RV_HARPS);
jitter     = jitter_proto;

% Linear part %
MJD    = MJD - min(MJD);
idx_linear  = (MJD < 0.024) | (MJD > 0.1);
MJD_linear = MJD(idx_linear);
RV_HARPS_linear    = RV_HARPS(idx_linear);
weight_linear       = weight(idx_linear);
fun_b = @(b) sum((b(1) + b(2) * MJD_linear - RV_HARPS_linear).^2 .* weight_linear);
b0 = [max(RV_HARPS), -max(RV_HARPS)/max(MJD)];
options = optimoptions(@fminunc,'Algorithm','quasi-newton');
b = fminunc(fun_b,b0,options);
           
            
h = figure;
    subplot(16,1,1:14)
    hold on 
    scatter(MJD*24, RV_FT, 30, 'rD', 'MarkerFaceColor', 'r', 'MarkerFaceAlpha', 0.5);
    scatter(MJD*24, RV_HARPS, 30, 'bo', 'MarkerFaceColor', 'b', 'MarkerFaceAlpha', 0.5);
%     errorbar(MJD2, RV_FT2, RV_noise2, 'k.', 'MarkerSize', 0.1);
%     errorbar(MJD2, RV_HARPS2, RV_noise2, 'b.', 'MarkerSize', 0.1);
    plot1 = plot(MJD*24, b(1) + b(2) * MJD, 'b--', 'LineWidth', 2);
    plot1.Color(4) = 0.4;        
    
    hold off
    xlim([-0.003*24 (max(MJD)+0.003)*24])
    legend('RV_{FT,L}', 'RV_{HARPS}', 'Linear trend')
%     title('Transit RV curve')
    ylabel('RV [m/s]')
    ylim([-59 49])
    set(gca,'fontsize', 15)
    set(gca,'xticklabel',[])    
    grid on

    % Jitter part %
    subplot(16,1,15:16)
    t2_plot         = linspace(min(MJD), max(MJD), 1001);
    jitter2_plot    = FUNCTION_GAUSSIAN_SMOOTHING(MJD, jitter, weight, t2_plot, width);
    hold on
    scatter(MJD*24, jitter - mean(jitter), 30, 'ks', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.5);
    errorbar(MJD*24, jitter - mean(jitter), RV_noise, 'k.', 'MarkerSize', 0.1)
    plot2 = plot(t2_plot*24, jitter2_plot - mean(jitter), 'k-', 'LineWidth',2);
    plot2.Color(4) = 0.4;       
    hold off
    xlim([-0.003*24 (max(MJD)+0.003)*24])
%     ylim([-3.2 3.2])
    xlabel('Time [h]')
    ylabel('\DeltaRV_{L} [m/s]')
    set(gca,'fontsize', 15)
    saveas(gcf,'8-Proto_jitter0','png')
close(h)

jitter_smooth  = FUNCTION_GAUSSIAN_SMOOTHING(MJD, jitter, weight, MJD, width);
RV_RM = RV_HARPS - (b(1) + b(2) * MJD);
fun_a = @(a) sum((a(1) * jitter_smooth + a(2) - RV_RM).^2 .* weight);
options = optimoptions(@fminunc,'Algorithm','quasi-newton');
a0 = [10, 0];
a = fminunc(fun_a,a0,options);

h = figure;
subplot(16,1,1:14)       % add first plot in 2 x 1 grid
    hold on 
    scatter(MJD*24, a(1) * jitter + a(2), 30, 'ks', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.5);
    plot3 = plot(t2_plot*24, a(1) * jitter2_plot + a(2) , 'k', 'LineWidth', 2);
    plot3.Color(4) = 0.4;  
    scatter(MJD*24, RV_RM, 30, 'bo', 'MarkerFaceColor', 'b', 'MarkerFaceAlpha', 0.5);    
%     jitter_final_Z = a(1) * jitter2 + a(2);
%     err_jitter_final_Z = RV_noise2 * sqrt(3) * a(1);
    errorbar(MJD*24, a(1) * jitter + a(2), (RV_noise * a(1)), 'k.', 'MarkerSize', 0.1)
    errorbar(MJD*24, RV_RM, RV_noise, 'b.', 'MarkerSize', 0.1)
    xlim([-0.003*24 (max(MJD)+0.003)*24])
    ylim([-85 85])
    ylabel('RV [m/s]')
    set(gca,'fontsize', 15)
    set(gca,'xticklabel',[])
    grid on 
    legend('\alpha \DeltaRV_{L}', 'WMA', 'RM effect')
    hold off
    
subplot(16,1,15:16)       % add first plot in 2 x 1 grid
    hold on
%     scatter(MJD2, a(1) * jitter2 + a(2) - RV_RM2, 30, 'ks', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.5)
%     errorbar(MJD2, a(1) * jitter2 + a(2) - RV_RM2, (RV_noise2 * sqrt(2) * a(1)), 'k.', 'MarkerSize', 0.1)
    scatter(MJD*24, a(1) * jitter_smooth + a(2) - RV_RM, 30, 'ko', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.5)
%     plot5 = plot(MJD2, a(1) * jitter_smooth2 + a(2) - RV_RM2, 'k', 'LineWidth', 2);
%     plot5.Color(4) = 0.4;
%     plot4 = plot(MJD2, a(1) * jitter_smooth2 + a(2) - RV_RM2 , 'k', 'LineWidth', 2);
%     plot4.Color(4) = 0.4;
    hold off
    xlabel('Time [h]')
    ylabel('Residual [m/s]')
%     ylim([-30 30])
    xlim([-0.003*24 (max(MJD)+0.003)*24])
    ylim([-12 12])
    set(gca,'fontsize', 15)
    saveas(gcf,'9-RM_fit0','png')
close(h)











%%%%%%%%%%%%%%%%%%
% Fitting Part 1 %   
%%%%%%%%%%%%%%%%%%
idx_t   = (MJD>54301.08) & (MJD<54301.241);
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
RV_HARPS2   = RV_HARPS2 - mean(RV_HARPS2);
% RV_FT2      = RV_FT(idx_t);
% RV_FT2      = RV_FT2 - mean(RV_FT2);
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
options = optimoptions(@fminunc,'Algorithm','quasi-newton');
b = fminunc(fun_b,b0,options);
           
            
h = figure;
    subplot(16,1,1:14)
    hold on 
%     scatter(MJD2*24, RV_FT2, 30, 'k+', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.5);
    scatter(MJD2*24, RV_HARPS2, 20, 'bo', 'MarkerFaceColor', 'b', 'MarkerFaceAlpha', 0.5);
%     errorbar(MJD2, RV_FT2, RV_noise2, 'k.', 'MarkerSize', 0.1);
%     errorbar(MJD2, RV_HARPS2, RV_noise2, 'b.', 'MarkerSize', 0.1);
    plot1 = plot(MJD2*24, b(1) + b(2) * MJD2, 'b--', 'LineWidth', 2);
    plot1.Color(4) = 0.4;        
    plotx = 0.025*24*ones(1,20);
    ploty = linspace(-59, 49, 20);
    plot(plotx, ploty, 'k--', 'Linewidth', 1)
    plotx = 0.106*24*ones(1,20);
    plot(plotx, ploty, 'k--', 'Linewidth', 1)    
    text(0.026*24, -10, 'ingress', 'FontSize', 15)
    text(0.107*24, -10, 'egress', 'FontSize', 15)    
    hold off
    xlim([-0.003*24 (max(MJD2)+0.003)*24])
    legend('RV_{HARPS}', 'Linear trend')
%     title('Transit RV curve')
    ylabel('RV [m/s]')
    ylim([-59 49])
    set(gca,'fontsize', 15)
    set(gca,'xticklabel',[])    
    grid on

    % Jitter part %
    subplot(16,1,15:16)
    t2_plot         = linspace(min(MJD2), max(MJD2), 1001);
    jitter2_plot    = FUNCTION_GAUSSIAN_SMOOTHING(MJD2, jitter2, weight2, t2_plot, width);
    hold on
    scatter(MJD2*24, jitter2 - mean(jitter2), 20, 'k^', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.5);
    errorbar(MJD2*24, jitter2 - mean(jitter2), RV_noise2*3^0.5, 'k.', 'MarkerSize', 0.1)
    plot2 = plot(t2_plot*24, jitter2_plot - mean(jitter2), 'r-', 'LineWidth',2);
    plot2.Color(4) = 0.4;       
    plotx = 0.025*24*ones(1,20);
    ploty = linspace(-10, 10, 20);
    plot(plotx, ploty, 'k--', 'Linewidth', 1)
    plotx = 0.106*24*ones(1,20);
    plot(plotx, ploty, 'k--', 'Linewidth', 1)       
    hold off
    xlim([-0.003*24 (max(MJD2)+0.003)*24])
%     ylim([-3.2 3.2])
    xlabel('Time [h]')
    ylabel('\DeltaRV_{H} [m/s]')
    set(gca,'fontsize', 15)
    grid on 
    saveas(gcf,'8-Proto_jitter2','png')
close(h)

jitter_smooth2  = FUNCTION_GAUSSIAN_SMOOTHING(MJD2, jitter2, weight2, MJD2, width);
RV_RM2 = RV_HARPS2 - (b(1) + b(2) * MJD2);
fun_a = @(a) sum((a(1) * jitter_smooth2 + a(2) - RV_RM2).^2 .* weight2);
options = optimoptions(@fminunc,'Algorithm','quasi-newton');
a0 = [10, 0];
a = fminunc(fun_a,a0,options);

h = figure;
subplot(16,1,1:14)       % add first plot in 2 x 1 grid
    hold on 
    scatter(MJD2*24, a(1) * jitter2 + a(2), 20, 'k^', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.5);
    plot3 = plot(t2_plot*24, a(1) * jitter2_plot + a(2) , 'r-', 'LineWidth', 2);
    plot3.Color(4) = 0.4;  
    scatter(MJD2*24, RV_RM2, 20, 'bo', 'MarkerFaceColor', 'b', 'MarkerFaceAlpha', 0.5);    
%     jitter_final_Z = a(1) * jitter2 + a(2);
%     err_jitter_final_Z = RV_noise2 * sqrt(3) * a(1);
    errorbar(MJD2*24, a(1) * jitter2 + a(2), (RV_noise2 * 3^0.5 * a(1)), 'k.', 'MarkerSize', 0.1)
%     errorbar(MJD2*24, RV_RM2, RV_noise2, 'b.', 'MarkerSize', 0.1)
    plotx = 0.025*24*ones(1,20);
    ploty = linspace(-70, 70, 20);
    plot(plotx, ploty, 'k--', 'Linewidth', 1)
    plotx = 0.106*24*ones(1,20);
    plot(plotx, ploty, 'k--', 'Linewidth', 1)   
    text(0.026*24, -60, 'ingress', 'FontSize', 15)
    text(0.107*24, -60, 'egress', 'FontSize', 15)
    xlim([-0.003*24 (max(MJD2)+0.003)*24])
    ylim([-70 70])
    ylabel('RV [m/s]')
    set(gca,'fontsize', 15)
    set(gca,'xticklabel',[])
    grid on 
    legend('Scaled \DeltaRV_{H}', 'WMA', 'RV_{HARPS} detrended')
    hold off
    
subplot(16,1,15:16)       % add first plot in 2 x 1 grid
    hold on
%     scatter(MJD2, a(1) * jitter2 + a(2) - RV_RM2, 20, 'ks', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.5)
%     errorbar(MJD2, a(1) * jitter2 + a(2) - RV_RM2, (RV_noise2 * sqrt(2) * a(1)), 'k.', 'MarkerSize', 0.1)
    scatter(MJD2*24, a(1) * jitter_smooth2 + a(2) - RV_RM2, 20, 'ko', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.5)
%     plot5 = plot(MJD2, a(1) * jitter_smooth2 + a(2) - RV_RM2, 'k', 'LineWidth', 2);
%     plot5.Color(4) = 0.4;
%     plot4 = plot(MJD2, a(1) * jitter_smooth2 + a(2) - RV_RM2 , 'k', 'LineWidth', 2);
%     plot4.Color(4) = 0.4;
    plotx = 0.025*24*ones(1,20);
    ploty = linspace(-12, 12, 20);
    plot(plotx, ploty, 'k--', 'Linewidth', 1)
    plotx = 0.106*24*ones(1,20);
    plot(plotx, ploty, 'k--', 'Linewidth', 1)   
    hold off
    xlabel('Time [h]')
    ylabel('Residual [m/s]')
%     ylim([-30 30])
    xlim([-0.003*24 (max(MJD2)+0.003)*24])
    ylim([-12 12])
    set(gca,'fontsize', 15)
    grid on 
    saveas(gcf,'9-RM_fit2','png')
close(h)



dlmwrite('MJD_2012.txt', MJD-MJD(1))
dlmwrite('Jitter_model_2012.txt', jitter_proto)
dlmwrite('RV_HARPS_2012.txt', (RV_HARPS - mean(RV_HARPS)) * 1000)
dlmwrite('RV_FT_2012.txt', (RV_FT - mean(RV_FT))' * 1000)










