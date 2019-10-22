% modified from 0816-FT-multiple_stars/FT.m

DIR         = '/Volumes/DataSSD/MATLAB_codes/0726-LSD_profile/';
file_list   = dir([DIR, 'LSDprof_norm/*.norm']);
file_name   = {file_list.name};
N_FILE      = size(file_name, 2);

% Obtain Julian dates
JDates      = importdata([DIR, 'JDates + RV shift']);
JD          = JDates.data(:,1);

% Obtain a first LSD line profiles
delimiterIn = ' ';
headerlinesIn = 2;
profD       = importdata([DIR, 'LSDprof_norm/', 'profD_lopeg_16aug14_v_02.s.norm'], delimiterIn, headerlinesIn);
x0          = profD.data(:,1) + 20;                                         % RV (20 km/s offset corrected)
A0          = profD.data(:,2);                                              
A0          = 1 - A0;
idx         = (x0>-150) & (x0<150);
x           = x0(idx);
A           = A0(idx);
grid_size   = mean(diff(x));                                                % 1.8
Fs          = 1/grid_size;

% estimate the size of array FFT_power
[aa, bb, yy]= FUNCTION_FFT(A, grid_size);
size1       = length(bb);
FFT_power   = zeros(size1, N_FILE);
Y           = zeros(size1, N_FILE);

%%%%%%%%%%%%%%%%%%%%%%%%
% Stacked LSD profiles %
%%%%%%%%%%%%%%%%%%%%%%%%
h = figure;
hold on
for n = 1:N_FILE

    profD       = importdata([DIR, 'LSDprof_norm/', char(file_name(n))], delimiterIn, headerlinesIn);
    A           = profD.data(:,2);
    A           = 1 - A(idx);
    plot(x, A)
    [FFT_frequency, FFT_power(:, n), Y(:, n)] = FUNCTION_FFT(A, Fs);
end     
hold off
% ylim([-0.001 0.1])
set(gca,'fontsize',20)
xlabel('km/s')
ylabel('Stockes I Intensity')
title('Normalized LSD profile')
saveas(gcf,'1-Normalized_LSD_profile','png')
close(h)

%%%%%%%%%%%%%%%%%%%%%%%%%%
% FT power in all epochs %
%%%%%%%%%%%%%%%%%%%%%%%%%%
h = figure; 
hold on
for n = 1:N_FILE
    plot(FFT_frequency, FFT_power(:, n), 'k')
%     title('Stacked power spectrum')
end 
hold off
xlabel('\xi [s/km]')
ylabel('Power')   
xlim([-0.03 0.03])
set(gca,'fontsize',20)
saveas(gcf,'2-FT_power','png')
close(h)

%%%%%%%%%%%%%%%
% Phase angle %
%%%%%%%%%%%%%%%
h = figure;
plot(FFT_frequency, unwrap(angle(Y(:, 1))), '.')
title('Phase angle')
xlabel('FT frequency (1 / velocity in wavelength)')
ylabel('Phase angle [radian]')
xlim([-0.35 0.35])
saveas(gcf,'3-Phase_angle','png')
close(h)

%%%%%%%%%%%%%%%%%%%%%
% Phase angle -> RV %
%%%%%%%%%%%%%%%%%%%%%
f_max = 0.008;
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
f_HL = 0.002;
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

RV_FT   = (RV_FT - mean(RV_FT))';
RV_FTL  = (RV_FTL - mean(RV_FTL))';
RV_FTH  = (RV_FTH - mean(RV_FTH))';
jitter_L  = RV_FT - RV_FTL;
jitter_H  = RV_FTH - RV_FT;

% RV_HARPS vs RV_FT %
h = figure; 
    hold on
    plot(RV_FT, jitter_L, 'b.', 'MarkerSize', 20)
    plot(RV_FT, jitter_H, 'r.', 'MarkerSize', 20)
    legend('\Delta RV_{L}', '\Delta RV_{H}', 'Location','northwest')
    hold off
    xlabel('RV_{FT} [km/s]')    
    ylabel('Jitter metrics [km/s]')
    saveas(gcf,'7-Jitter','png')
close(h)

% Time series 
figure; 
    plot(JD, RV_FT*1000,  'r.', 'MarkerSize', 20)
    xlabel('Time [d]')
    ylabel('RV [m/s]')

dlmwrite('RV_FT.txt', RV_FT)
dlmwrite('RV_FTL.txt', RV_FTL)
dlmwrite('RV_FTH.txt', RV_FTH)

% [DIR, 'LSDprof_norm/', char(file_name(1))]