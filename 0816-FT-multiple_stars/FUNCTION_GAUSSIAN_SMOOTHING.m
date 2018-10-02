function [y_smooth] = FUNCTION_GAUSSIAN_SMOOTHING(t, y, weight2, t_smooth, len_smooth)

y_smooth = zeros(size(t_smooth));

for i = 1:length(t_smooth)
    weight = exp(-(t_smooth(i)-t).^2/(2*len_smooth^2)) .* weight2;
    y_smooth(i) = sum(y .* weight) / sum(weight);
end
