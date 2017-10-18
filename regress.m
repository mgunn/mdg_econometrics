function [results, options] = regress(varargin),
% Author: Matthew Gunn
% Notes:  In general, behavior is designed to match that of Stata
%
%                               ** USAGE ** 
% [results, options_used] = regress(y, X, <option>, <option_val>, ...)
%
% options:                   values:
%          'vce'                  can take 'ols', 'robust', or 'cluster'    (default 'ols')
%          'cluster_var'          n by 1 vector or string cell array of ids (default blank)
%          'noconstant'           either 1 or 0                             (default 0)
%
% Example: [r, o] = mdg_econometrics.regress(d.math4, [d.lunch, d.expp], 'vce','cluster',
%                                            'cluster_var', d.distid, 'noconstant', 1);
%
% Output: results.b        the constant term will be the LAST entry
%         results.bse      standard errors
%         results.bcov     covariance matrix
%         results.rss      residual sum of squares
%         results.mss      model sum of squares
%         results.r2   
%         results.F 
%         results.Fval     pulled from F cdf
%         results.t        t-stat
%         results.df       degrees of freedom (n-k)
%         results.pvals    pvals for results.b
%         results.conf95   95% confidence interval

%%                         PREP WORK                             %%

options = struct();
options.command = 'regress';

% set defaults
y = varargin{1};
X = varargin{2};

if(islogical(X))
    X = double(X);
end

if(~isempty(X))
    if(size(y,1) ~= size(X,1))
        error('rows of y and X are different');
    end
end

options.vce = 'ols';
options.noconstant = false;



%%                   Set user defined options                   %%
if(nargin < 2)
    error('[results, options_used] = regress(y, X, ...)');
end

if(mod(nargin, 2) ~= 0)
    error('Should be even number of arguments');
end

i = 3;
while(i <= nargin)
    switch varargin{i}
      case 'vce'
        options.vce = varargin{i+1};
      case 'cluster_var'
        options.cluster_var = varargin{i+1};
      case 'noconstant'
        options.noconstant = varargin{i+1};
        if(~islogical(options.noconstant))
            error('argument for noconstant option should be logical: i.e. true or false');
        end
      otherwise
        error(sprintf('Option ''%s'' not understood.',varargin{i}));
    end
    i = i + 2;
end

if(~options.noconstant)
    X = [X, ones(length(y), 1)];
end

%%                      SETUP                 %%
if(size(y, 2) > 1)
  error('only supports 1 dimensional y');
end;

if(isfield(options,'cluster_var'))
    mask_col = ~any(isnan([X, y, options.cluster_var]),2);
else
    mask_col = ~any(isnan([X, y]),2);
end

X = X(mask_col,:);
y = y(mask_col,:);
if(isfield(options,'cluster_var'))
    options.cluster_var = options.cluster_var(mask_col, :);
end


%%                      Additional checks                       %%
if(strcmp(options.vce,'cluster'))
    if(~isfield(options,'cluster_var'))
        error('If vce=''cluster'' is set, then ''cluster_var'' must be set to an n by 1 vector or cell array of strings');
    else
        [n3, shouldbe1] = size(options.cluster_var);
        if(n3 ~= length(y))
            error('Number of rows of cluster_var should be same as y and X');
        end
        if(shouldbe1 ~= 1)
            error('Only 1 cluster id');
        end        
    end
end

% $$$ temp = mdg_econometrics.removenan([y, X]);        %% REMOVE NaN rows!!!!
% $$$ y = temp(:,1);
% $$$ X = temp(:,2:end);
% $$$ 

[n, cols] = size(y);
[n2, k] = size(X);


%%                         ESTIMATION                            %%
XprimX     = X'*X;
Xprimy     = X'*y;


% Formula 1:
% results.b  = X \ Y;           % solve b = (x'x)^-1 (x' y)

% Formula 2:     reuse XprimX for later
% results.b     = XprimX \ Xprimy;        

% Formula 3:     reuse cholesky decomposition for later so we can just backsolve (Chol and backsolve is O(n^2) op)
R = chol(XprimX);   % Cholesky decomposition:    XprimX_chol_UT' * XprimX_chol_UT = XprimX
                    % 
                    % We want to solve: XprimX * b = Xprimy
                    %                     R'*R * b = Xprimy
                    %                            b = R\(R'\Xprimy)

results.b = R\(R'\Xprimy);


%%                     STANDARD ERRORS                          %%
yhat  = X * results.b;

u     = y - yhat;   %calculate residual
results.resid = u;               %store for later

%% Deriviation of formula %%
%  var(B) = var(inv(X'X)(X' (X*BETA + u)))
%         = var( inv(X'X) (X'u))
%         = inv(X'X) var(X'U) inv(X'X) 
%         = [inv(X'X/n) ] var(X'U/n) [inv(X'X/n)]
% X'X/n         --p-->  M
% (n^-.5) X'U   --d-->  N(0, V)
% V=E[u^2 x*x']
%
%  =>var(X'U/n)=var(n^-.5 X'U)/n=V/n
% 
%
%  var(B) = BETA + M^-1 V/n M^-1

%% finish calculation of asymptotic covariance matrix

switch options.vce
  case 'ols'
    s2 = u' * u / (n - k);
    results.bcov = inv(XprimX) * s2;
  
  case 'robust'
    temp = (u * ones(1, k)) .* X;
    V    = temp' * temp / n;
    
    % Recall R'*R = XprimX by cholesky decomposition
    % COV(b) = M^-1 * V * M^-1
    %        = n^2 * XprimX^-1 * V * XprimX^ - 1
    %        = n^2 * (R'*R)^-1 * V * (R'*R)^-1
    %        = n^2 * R^-1 * R'^-1 * V * R^-1 * R'^-1
    %    results.bcov = inv(XprimX) * V * inv(XprimX) * (n^2 / (n - k));
    results.bcov = (R \ (R' \ ((V / R) / R'))) * (n^2 / (n - k));           %reuses Cholesky decomposition
  case 'cluster'    
    [V, n_c] = mdg_econometrics.getClusteredSMatrix(u, options.cluster_var, X);
    results.bcov = (R \ (R' \ ((V / R) / R'))) * (n * ((n - 1) / (n - k)) * (n_c / (n_c - 1))); 
    % Notes:
    % The n is just to deal with that V is normalized by n
    %
    % The ((n - 1) / (n - k)) * (n_c / (n_c - 1))) is adjustment. Used by Stata, Hansen code etc...
    results.n_clusters = n_c;

  case 'cluster2'
    if(size(options.cluster_var,2)==2)
        cluster_id1 = options.cluster_var(:,1);
        cluster_id2 = options.cluster_var(:,2);
        [~, ~, cluster_id12] = unique(options.cluster_var,'rows');        
    else
        error('clusterids has bad size');
    end
       
    [S1, n_c1] = mdg_econometrics.getClusteredSMatrix(u, cluster_id1, X);  %% edit 2-25-2017 taking out w, is this right?
    [S2, n_c2] = mdg_econometrics.getClusteredSMatrix(u, cluster_id2, X);
    [S12, n_c12] = mdg_econometrics.getClusteredSMatrix(u, cluster_id12, X);

% $$$     [S1, n_c1] = mdg_econometrics.getClusteredSMatrix(u, cluster_id1, Xw);  %% edit 2-25-2017 taking out w, is this right?
% $$$     [S2, n_c2] = mdg_econometrics.getClusteredSMatrix(u, cluster_id2, Xw);
% $$$     [S12, n_c12] = mdg_econometrics.getClusteredSMatrix(u, cluster_id12, Xw);

    
    
% $$$     V1  = (R \ (R' \ ((S1 / R) / R'))) * (n * ((n - 1) / (n - k)) * (n_c1 / (n_c1 - 1)));
% $$$     V2  = (R \ (R' \ ((S2 / R) / R'))) * (n * ((n - 1) / (n - k)) * (n_c2 / (n_c2 - 1)));
% $$$     V12 = (R \ (R' \ ((S12 / R) / R'))) * (n * ((n - 1) / (n - k)) * (n_c12 / (n_c12 - 1)));
    V1  = (R \ (R' \ ((S1 / R) / R'))) * n;
    V2  = (R \ (R' \ ((S2 / R) / R'))) * n;
    V12 = (R \ (R' \ ((S12 / R) / R'))) * n;
    
    V = V1 + V2 - V12;
    results.bcov = V;
    %    results.bcov =  
    % Notes:
    % The n is just to deal with that V is normalized by n
    %
    % The ((n - 1) / (n - k)) * (n_c / (n_c - 1))) is adjustment. Used by Stata, Hansen code etc...
    %    results.n_clusters = n_c;    
    results.n_c1 = n_c1;
    results.n_c2 = n_c2;
    results.n_c12 = n_c12;
    
    
  otherwise
    error('Unknown VCE but code should have already checked for this!?!?');
end

results.rss  = u'*u;                                   % residual sum of squares
results.mss  = norm(yhat - mean(y), 2).^2;             % model sum of squares

results.r2    = results.mss / (results.rss + results.mss);                           %n / var(y,1);


if(strcmp(options.vce, 'cluster'))
    results.df_r = n_c - 1;
elseif(strcmp(options.vce, 'cluster2'))
    results.df_r = NaN; % FIXME: What should go here?
else
    results.df_r = n - k;
end
if(~options.noconstant) % if has a constant
    results.df_m = k - 1;
    if( k - 1 == 0)
        results.F = NaN;
    else
        results.F     = results.b(1:end-1)' * inv(results.bcov(1:end-1, 1:end-1)) * results.b(1:end-1) / (k-1);
    end
else                   % if does not have a constant
    results.df_m = k;
    results.F     = results.b' * inv(results.bcov) * results.b / k;
end
results.Fval = 1 - fcdf(results.F, results.df_m, results.df_r);

results.bse = diag(results.bcov).^.5;

results.t = results.b ./ results.bse;


%%%% FIXME: df_r should be n - k always.... this df_m has it munged up for F test

if(strcmp(options.vce, 'cluster2'))
    results.pvals = NaN(size(results.b));
    results.pvals(results.t >= 0) = 2 * (1 - normcdf(results.t(results.t>=0), 0, 1));
    results.pvals(results.t < 0)  = 2 * normcdf(results.t(results.t<0), 0, 1);
else
    results.pvals = NaN(size(results.b));
    results.pvals(results.t >= 0) = 2 * (1 - tcdf(results.t(results.t>=0), results.df_r));
    results.pvals(results.t < 0)  = 2 * tcdf(results.t(results.t<0), results.df_r);
end
    
%results.R2adj = 1 - results.s2 ./ var(y);

results.conf95 = [results.b + results.bse * tinv(.025, results.df_r), results.b + results.bse * tinv(.975, results.df_r)];
results.n_obs  = n;