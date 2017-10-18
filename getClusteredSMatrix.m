function [V, n_c] = getClusteredSMatrix(resid, clusterid, X)
%
% This function is used behind the scenes by ivregress and others
%
% Requires: mdgtools.mg_getRowsWithKey
%
% [V, n_c] = getClusteredSMatrix(resid, clusterid, X)
%
%              resid: vector of residuals 
%          clusterid: vector of cluster ids or cell array of unique strings identifying cluster
%                  X: right hand side matrix from regression
%
%                  V: the matrix we're estimating 
%                n_c: number of clusters

%%                                          START THE STUFF                                            %%
% c_unique = unique(c);                                     % METHOD 1: slow as shit and wasteful; use METHOD 2

[m, n] = size(clusterid);
if(n~=1)
    error('clusterid has more than 1 column...');
end
if(m ~= size(X,1))
    error('rows of clusterid and X dont match');
end

[c_unique, keyrowmap] = mdgtools.mg_getRowsWithKey(clusterid);      % METHOD 2: MUCH faster
                                                                    %  c++ function that calculates this FAST

n_c = length(c_unique);

[n_obs, k] = size(X);

V = zeros(k,k);    
for i=1:n_c
    %        rows = (c==c_unique(i));         % Method 1: use a mask (requires above method 1 code uncommented)
    rows = keyrowmap{i};                      % Method 2: NOT a mask. Just a list of rows
    
    X_ci = X(rows,:);                         % Select relevant section of X matrix
    
    resid_ci = resid(rows);                   % select relevant section of resids vector
    n_ci = length(resid_ci);                  % elements of cluster i
    
    %    temp = X_ci .* (resid_ci * ones(1,k));
    %    V_ci = temp' * temp / n;

    Xprim_resid_ci = X_ci' * resid_ci;
    
    V = V + Xprim_resid_ci * Xprim_resid_ci';
    
    %    V_ci = X_ci' * resid_ci * resid_ci' * X_ci;
    %    V = V + V_ci;
    
    %    M = M + M_ci;
end

V = V / n_obs;