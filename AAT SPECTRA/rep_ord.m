function [wav_arr, int_arr, err_arr, int_nor_arr] = rep_ord(A, wav_min, wav_max)
% returns the array (interpolated) -> replace the spectrum order

% set up region for order 19
A1 = A(7:end-1);
M  = double(zeros(size(A1,2), 3));

for i = 1 : size(A1,2)
    M(i,:) = str2num(A1{i});     % Kein Problem mit der Warnung. 
end

idx = or((M(:,3) == 33), (M(:,3) == 31));
M = M(idx,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test if there's outlier, e.g. cosmic rays %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wav         = M(:,1);
int         = M(:,2);
idx         = (wav > wav_min) & (wav < wav_max);
wav         = wav(idx);
int         = int(idx);

DEG_CLEAN   = 1;
DEG_POLY    = 1;
[X, Y, ~]   = rep_nor(wav, int, DEG_CLEAN, DEG_POLY); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normalization after getting rid of outliers %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 2 works slightly better than 3. 
% 1 works better than 2 for shorter wavelength range
DEG_CLEAN   = 1;                
% 2 performs better than 1, however 1 is more robust, especially in the
% precense of cosmic rays.
DEG_POLY 	= 1;                
[wav, int, int_nor_2]  = rep_nor(X, Y, DEG_CLEAN, DEG_POLY);

%%%%%%%%%%%%%%%%%%%%
% 1D interpolation %
%%%%%%%%%%%%%%%%%%%%
wav_arr     = (wav_min : 0.01 : wav_max)';
int_arr     = interp1q(wav, int, wav_arr);
err_arr     = int_arr.^0.5;
int_nor_arr = interp1q(wav, int_nor_2, wav_arr);