function [X, Y, Y_nor] = rep_nor(wav, int, DEG_CLEAN, DEG_POLY)

x   = wav;
y   = int;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First cleanup - get rid of large fluctuations %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
max_diff = max(diff(y));

for i = 1:10
    idx = [(abs(diff(y)) < max_diff * 0.1)', 0]';        % index to be kept
    idx = idx | circshift(idx, 1);      % circshift move the array rightward or downward
    if sum(1 - idx) == 0
        break
    end
    x   = x(idx);
    y   = y(idx);
end

%%%%%%%%%%%%%%%%%%%%%%%
% locate local maxima %
%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:DEG_CLEAN
    [x, y] = local_max(x, y);
end

%%%%%%%%%%%%%%%%%%%%%%%
% remove local minima %
%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:20
    [x_min, y_min] = local_min(x, y);
    if isempty(x_min) && isempty(y_min)
        break
    end
    idx = (~ismember(x, x_min)) & (~ismember(y, y_min));
    x   = x(idx);
    y   = y(idx);
end

%%%%%%%%%%%%%%%%%%%
% Standardization %
%%%%%%%%%%%%%%%%%%%
MEAN    = mean(x);
STD     = std(x) ;
x       = (x - MEAN) / STD;

%%%%%%%%%%%%%%%%%%
% Polynomial fit %
%%%%%%%%%%%%%%%%%%
P1  = polyfit(x, y, DEG_POLY);
y_fit   = polyval( P1, (wav-MEAN)/STD );
y_nor   = int ./ y_fit;

%%%%%%%%%%%%%%%%%%%%%%%
% Get rid of outliers %
%%%%%%%%%%%%%%%%%%%%%%%
% Outliers correction is still required at this stage. 
torr        = 0.1;
idx_torr    = y_nor<(1+torr);
X           = wav(idx_torr);
Y           = int(idx_torr);
Y_nor       = y_nor(idx_torr);
