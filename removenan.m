function [result mask_col] = removenan(matrix);
%  removenan(matrix)
%  
%  Author: Matthew Gunn
%       Remove rows containing NaN
%
%       result is the resulting matrix
%      
%       mask_col is a vector with the same number of rows as matrix
%       such that every entry of mask_col is TRUE if row is
%       retained and FALSE if the row is removed.

mask_col = ~any(isnan(matrix),2);
result = real(matrix(mask_col,:)); %edit 7/4/2017 added call to real

% $$$ [rows, cols] = size(matrix);
% $$$ 
% $$$ % $$$ I don't understand why second col?
% $$$ % $$$ temp = isnan(matrix) * ones(cols, cols);
% $$$ % $$$ mask_col = temp(:,1)==0;
% $$$ 
% $$$ temp = isnan(matrix) * ones(cols, 1);
% $$$ mask_col = temp==0;
% $$$ 
% $$$ 
% $$$ mask = (mask_col * ones(1, cols)) == 1;
% $$$ T = sum(mask_col);
% $$$ 
% $$$ %%%%%%% pull it out!
% $$$ result = reshape(matrix(mask), T, cols);
