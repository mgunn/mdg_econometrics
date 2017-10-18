% function output_data = winsorize(input_data, level)
% input:
% level:  eg. .01 or .05
function output_data = winsorize(input_data, level)
if(level > .2)
    error('really, level above .2?');
end
new_max = quantile(input_data, 1 - level);
new_min = quantile(input_data, level);

k = size(input_data,2);
output_data = input_data;

for i=1:k
    output_data(input_data(:,i) < new_min(i), i) = new_min(i);
    output_data(input_data(:,i) > new_max(i), i) = new_max(i);
end