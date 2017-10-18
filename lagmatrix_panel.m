% function X_out = lagmatrix_panel(X_in, lags, ids, timevar)
%
% Author: Matthew Gunn
% 
% lagmatrix for panel data! This applies lagmatrix but groups by ids
function X_out = lagmatrix_panel(X_in, lags,ids)

[key_unique, keymap] = mdgtools.mg_getRowsWithKey(ids);

X_out = zeros(size(X_in));

for i=1:length(key_unique)
    X_out(keymap{i},:) = lagmatrix(X_in(keymap{i},:),lags);    
end
