% Modify from Project180131-FT_SOAP in order to pick up radnom phases
% instead of consecutive observations

% Use simulated spectra with planets AND stellar jitter: 
% /Volumes/DataSSD/SOAP_2/outputs/02.01/

%%%%%%%%%%%%%%
% Parameters %
%%%%%%%%%%%%%%
SN              = 5000;
N_FILE          = 100;                               
grid_size       = 0.1;
Fs              = 1/grid_size;
v0              = (-20 : grid_size : 20)';          % km/s
dir1            = '/Volumes/DataSSD/SOAP_2/outputs/02.01/';
dir2            = '/Volumes/DataSSD/SOAP_2/outputs/02.01/CCF_dat/';
RV              = importdata([dir1, 'RV.dat']) / 1000;      % activity induced RV [km/s]
idx             = (v0 > -10) & (v0 < 10);
v1              = v0(idx);
period_min      = 2;
period_max      = 100;

%%%%%%%%%%%%%%%%%%%
% Calculate Power %
%%%%%%%%%%%%%%%%%%%

% estimate the size of array FFT_power
filename    = [dir2, 'CCF', num2str(1), '.dat'];
A           = 1 - importdata(filename);
A           = A + normrnd(0, (1-A).^0.5/SN);
A           = A(idx);
A1          = A;
[aa, bb,yy]  = FUNCTION_FFT(A, Fs);
size1       = length(bb);
FFT_power   = zeros(size1, N_FILE);
Y           = zeros(size1, N_FILE);

% randomly select 200 different integers from 0 to 999
% 100 corresponds to one solar roation period ~ 25 days 
MJD     = sort(randperm(1000,N_FILE))-1;
jitter  = zeros(1,N_FILE);

jitter_full = zeros(1,1000);
for n = 0:999
    jitter_full(n+1) = RV(mod(n,100)+1);
end

RV_noise = zeros(1,N_FILE);

figure;
hold on

for n = 1:N_FILE
        
    i           = mod(MJD(n), 100);
    jitter(n)   = RV(i+1);
    v_planet    = 0 * sin(i/100*7*2*pi + 1) * 0.001;         % km/s
    filename    = [dir2, 'CCF', num2str(i), '.dat'];
    A           = 1 - importdata(filename);
    plot(v0(idx),A(idx)-A1, '.')
    A_spline    = spline(v0, A, v1+v_planet);
    A_spline    = A_spline + normrnd(0, (1-A_spline).^0.5/SN);
    f_tpl       = fit( v1, A_spline, 'a*exp(-((x-b)/c)^2)+d', 'StartPoint', [0.5, 0, 4, 0] );
    RV_noise(n) = f_tpl.b;

    [FFT_frequency, FFT_power(:, n), Y(:, n)] = FUNCTION_FFT(A_spline, Fs);

    if 1            % Plot line profile in frequency space
        h_fft = figure;
        hold on 
        plot(FFT_frequency(1:int16(size1/3)), FFT_power(1:int16(size1/3), n))
        plot(FFT_frequency(1:int16(size1/3)), FFT_power(1:int16(size1/3), n), '.', 'markers',12)
        hold off
        xlabel('FT frequency domain')
        ylabel('Normalized Power')        
        title_name = ['FT - file', num2str(n)];
        out_eps = [title_name, '.eps'];
        print(out_eps, '-depsc')
        close(h_fft); 
    end
    
end     
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct a template in FT space %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FFT_tpl = mean(FFT_power');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Correct each observation %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Y_c = zeros(size(Y));
Y_i = zeros(size(Y));
iFFT_power = zeros(size(FFT_power));

% Scale the power, and keep the phase

for n = 1:N_FILE
    for i = 1:length(FFT_tpl)
        Y_c(i,n) = Y(i,n) / abs(Y(i,n)) *  FFT_tpl(i);
    end
end

for n = 1:N_FILE
    [iFFT_power(:, n), Y_i(:,n)] = FUNCTION_iFFT(Y_c(:,n));
end

figure;
hold on
for n = 1:N_FILE
    plot(FFT_frequency, iFFT_power(:, n) - iFFT_power(:, 1), '.')
end
hold off 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Examine the corrected line centre %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_planet = zeros(1,N_FILE);
for n = 1:N_FILE
    % idx     = iFFT_power(:, n) > 0.01;
    % f_power = fit( FFT_frequency(idx)'-0.5, iFFT_power(idx, n), 'a*exp(-((x-b)/c)^2)+d', 'StartPoint', [1, 0.5, 1, 0] );
    idx = (FFT_frequency<1) & (FFT_frequency>0);
    f_power = fit( FFT_frequency(idx)', iFFT_power(idx, n), 'a*exp(-((x-b)/c)^2)+d', 'StartPoint', [0.5, 0.5, 0.18, 0] );
    i_planet(n) = f_power.b;
end

% test %
figure; plot(f_power, FFT_frequency(:), iFFT_power(:, n))

v_planet    = 2 * sin(mod(MJD, 100)/100*7*2*pi + 1) * 0.001;
figure; plot(RV_noise, i_planet, '.', 'markers',12)

