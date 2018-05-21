% Convert a Gaussian pulse from the time domain to the frequency domain.

%%%%%%%%%%%%%%%%%
% Generate data % 
%%%%%%%%%%%%%%%%%
Fs  = 100;              % Sampling frequency
t   = -0.5:1/Fs:0.5;    % Time vector 
L   = length(t);        % Signal length
% X   = 1/(4*sqrt(2*pi*0.01))*(exp(-t.^2/(2*0.01)));
X   = 1/(4*sqrt(2*pi*0.01))*(exp(-t.^2/(2*0.01)) + exp(-(t-0.1).^2/(1*0.01)));

if 1                    % Plot the pulse in the time domain
    f1 = figure;
    plot(t,X)
    title('Gaussian Pulse in Time Domain')
    xlabel('Time (t)')
    ylabel('X(t)')
end

%%%%%%
% FT % 
%%%%%%
[fshift,power] = FUNCTION_FFT_iFFT(X, Fs);

h = figure;
hold on 
plot(fshift,power) 
plot(fshift,power, '.', 'markers',12) 
hold off
title('Gaussian Pulse in Frequency Domain')
xlabel('Frequency (f)')
ylabel('|P(f)|')


%%%%%%%%
% FTFT % 
%%%%%%%%
% Fs  = 64;
% t   =





[fshift,power] = FUNCTION_FFT_noise(X, Fs);

f2 = figure;
hold on 
plot(fshift,power) 
plot(fshift,power, '.', 'markers',12) 
hold off
title('Gaussian Pulse in Frequency Domain')
xlabel('Frequency (f)')
ylabel('|P(f)|')