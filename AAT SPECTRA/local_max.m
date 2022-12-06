function [max_x,max_y] = local_max(x,y)

if length(x) ~= length(y)
    error('Two input data should have the same length.');
end

if (nargin < 2)||(nargin > 3)
    error('Please see help for INPUT DATA.');
end

% Find the extreme maxim values 
% and the corresponding indexes
%----------------------------------------------------
extrMaxIndex =   find(diff(sign(diff(y)))==-2)+1;
max_x = x(extrMaxIndex);
max_y = y(extrMaxIndex);