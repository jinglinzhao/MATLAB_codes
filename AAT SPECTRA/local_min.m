function [min_x,min_y] = local_min(x,y)

if length(x) ~= length(y)
    error('Two input data should have the same length.');
end

if (nargin < 2)||(nargin > 3)
    error('Please see help for INPUT DATA.');
end
    
% Find the extreme minim values 
% and the corresponding indexes
%----------------------------------------------------
extrMinIndex =   find(diff(sign(diff(y)))==+2)+1;
min_x = x(extrMinIndex);
min_y = y(extrMinIndex);