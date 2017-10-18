function [results, options] = xtreg(varargin),
% Author: Matthew Gunn
% Notes:  In general, behavior is designed to match that of Stata
%
%                               ** USAGE ** 
% [results, options_used] = xtreg(y, X, category_ids, <option>, <option_val>, ...)
%
% Results is a struct with results. options_used is a struct showing what xtreg actually used
%
% options:                   values:
%          'fe' (default)
%          'ra'            
%
%
%          'vce'                  can take 'ols', 'robust', or 'cluster'    (default 'ols')
%          'cluster_var'          n by 1 vector or string cell array of ids (default blank)
%          'noconstant'           either 1 or 0                             (default 0)
%
% Example: [r, o] = mdg_econometrics.xtreg(d.math4, [d.lunch, d.expp], d.distid, 'vce','cluster',
%                                            'cluster_var', d.distid);
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
%         results.df_r     degrees of freedom for residuals
%         results.pvals    pvals for results.b
%         results.conf95   95% confidence interval


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                           PROCESS OPTIONS and SETUP                             %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y = varargin{1};
X = varargin{2};
category_ids = varargin{3};

%warning('FIXME category id stuff and NAN is broken');

% Explanation of behavior:
% 

if(iscell(category_ids))    % idea: if it's a string, we don't want to check isnan
    index_nonan = ~any(isnan([y, X]), 2);
    y = y(index_nonan,:);
    X = X(index_nonan,:);
    category_ids = category_ids(index_nonan,:);    
else                       % if NOT a string, we want to check for nan
    index_nonan = ~any(isnan([y, X, category_ids]), 2);
    y = y(index_nonan,:);
    X = X(index_nonan,:);
    category_ids = category_ids(index_nonan,:);    
end

results = struct();

if(nargin < 3)
    error('[results, options_used] = xtreg(y, X, category_ids, <option>, <optino_val>...)');
end

% $$$ if(mod(nargin+1, 2) ~= 0)
% $$$     error('Should be odd number of arguments');
% $$$ end

options = struct();
options.command = 'xtreg';
options.estimator = 'fe';
options.vce = 'ols';

i = 4;
while(i <= nargin)
    switch varargin{i}
      case 'fe'
        options.estimator = 'fe';
        i = i + 1;
      case 'be'
        options.estimator = 'be';
        i = i + 1;
      case 'ra'
        options.estimator = 'ra';
        i = i + 1;
      case 'vce'
        options.vce = varargin{i+1};
        i = i + 2;
      case 'cluster_var'
        options.cluster_var = varargin{i+1};
        options.cluster_var = options.cluster_var(index_nonan);
        i = i + 2;
      case 'noconstant'
        options.noconstant = varargin{i+1};
        i = i + 2;
      otherwise
        error(sprintf('Option ''%s'' not understood.',varargin{i}));
    end
end

%%                             SETUP                            %%
[n, cols] = size(y);
if(cols > 1)
  error('only supports 1 dimensional y');
end;

[n2, k] = size(X);

if(n2 ~= n)
    error('rows of y and X are different');
end

%%                      Additional checks                       %%


if(strcmp(options.vce,'cluster'))
    if(~isfield(options,'cluster_var'))
        error('If vce=''cluster'' is set, then ''cluster_var'' must be set to an n by 1 vector or cell array of strings');
    else
        [n3, shouldbe1] = size(options.cluster_var);
        if(n3 ~= n)
            error('Number of rows of cluster_var should be same as y and X');
        end
        if(shouldbe1 ~= 1)
            error('Only 1 cluster id');
        end        
    end
end

%category_unique = unique(category);
%new_ids = mdgtools.getNewIds(category);

[n2, k] = size(X);

% Z = [X, D]
% [b; a] = (Z'Z)^-1(Z'Y)
% [b; a] = ([X'; D'] [X, D])^-1 [X'Y D'Y]
%        = (X'X X'D; D'X D'D)^-1 [X'Y D'Y]
%    etc....

if(n ~= length(category_ids))
    error('category ids have wrong length?');
end

% $$$ if(k >= 1 && norm(X(:,1) - ones(n,1)) < .0000000001)
% $$$     error('Looks like you included a constant term.  Thats wrong!');
% $$$ end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                               ESTIMATION                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First step is basically to demean vector y and matrix X based on category id
% 
% Md is a demeaning matrix 
% MdY is demeaned y vector
% MdX is demeaned X matrix


%% See Green p. 287
%formulas
%Md = eye(K) - D * inv(D' * D) *D';
%b  = inv(X' * Md * X) * (X' * Md * y);
%a  = inv(D' * D) *D' *(y - X*b)

% $$$ %n_cat = max(sequential_category_ids);
% $$$ %D = mdgtools.makeIndicators(sequential_category_ids);
% $$$ % method 1 to calcualte Md
% $$$ temp = (D' * D) \ D';         %solve inv(D' * D) * D'
% $$$ Md = eye(n) - D * temp;

%%                      CONSTRUCT WITHIN AND BETWEEN DATA                       %%
% method 2: (idea: giant matrices are slower than using loops sensibly)
[uniqueid, map_seqid_to_catid] = mdgtools.mg_getRowsWithKey(category_ids);
n_cat = length(uniqueid);

%% within data  %%
MdX = zeros(n, k);                 % MdX: Whole X matrix with mean subtracted for each category_id
Mdy = zeros(n, 1);                 % Mdy: Whole y vector with mean subtracted for each category_id

%% between data %%
ybar = zeros(n_cat, 1);            % Vector where each entry is the mean y for the category
Xbar = zeros(n_cat, k);            % Matrix where each entry is the mean X(i,:) for the category


for i=1:length(uniqueid)
    list = map_seqid_to_catid{i};
    n_stu_i = length(list);

    Xtemp = X(list,:);
    Xbar_i = mean(Xtemp, 1);
    MdX(list,:) = Xtemp - ones(n_stu_i,1) * Xbar_i;
    
    ytemp = y(list,:);
    ybar_i = mean(ytemp, 1);
    Mdy(list) = ytemp - ybar_i;
    
    ybar(i)   = ybar_i;
    Xbar(i,:) = Xbar_i;
end

switch options.estimator
    
  case 'fe'
    % (X' Md X)^-1 * (X' Md y)    
    % Algorithm 1:  ** SLOW **
    %Xprim_Md = X' * Md;
    %results.b = (Xprim_Md * X) \ (Xprim_Md * y);
    
    % Algorith 2:   ** FAST **
    % results.b = (MdX' * X) \ (MdX' * y);
    
    % Algorithm 3:
    XprimMdX = X' * MdX;
    XprimMdy = MdX' * y;           
    R        = chol(XprimMdX);   % Cholesky factorization allows us to quickly solve any linear system
                                 % involving XprimMdX whenever we want
    
    results.b = R\(R'\XprimMdy); % Use Cholesky decomposition to do two backsolves
    
    e         = Mdy - MdX * results.b;      % These are the residuals
    u         = y - Mdy;                    % These are the fixed effects
    ubar      = ybar - Xbar*results.b;
    
    %Note:
    % R'R = A                Cholesky returns R such that left statement is true. Then:
    % (R')^-1 * R^-1 = A^-1
    % x = A^-1 * b 
    %   = (R')^-1 * R^-1 * b
    %   = (R' \ (R \ b)
    
  case 'be'          %% FIXME: THIS IS TOTAL HACK %% 
    warning('Fstats and a bunch of things are wrong!!!!  Right now, between effects code is a total hack');
    
    %% FIXME: TOTAL HACK.  I'm renaming crap everywhere to do what I want later...
    
    Xbar  = [Xbar, ones(length(Xbar), 1)];
    MdX   = [X, ones(n, 1)];
    X     = [X, ones(n, 1)];
    R = chol(Xbar'*Xbar);

    results.b = R\(R'\(Xbar' * ybar))    

    
    %    results.b
    ubar = ybar - Xbar * results.b;
    e = ubar;

    %    results.b = results.b(1:end-1);
    
    n = length(Xbar); n_cat = 1;
    
    XprimMdX = Xbar' * Xbar;
    

    
% $$$     u = zeros(n, 1);
% $$$     
% $$$     for i=1:length(uniqueid)
% $$$         list = map_seqid_to_catid{i};
% $$$         u(list) = ubar(i);
% $$$     end
% $$$ 
% $$$     e = y - X * results.b - u;
    
    
    %    error('estimator ''be'' not supported yet');
    
  case 'ra'
    error('estimator ''ra'' not supported yet');
    
  otherwise
    error(sprintf('unknown estimator ''%s''',options.estimator));
end
    
%b  = x \ y;           %solve b = (x'x)^-1 (x' y)
%resid = y - x * b;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                              STANDARD ERRORS                                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deriviation of formula %%
%  var(b) = var(inv(x'x)(x' x*beta + e))
%         = beta + var( inv(x'x) (x'e))
%         = beta + inv(x'x) var(x'e) inv(x'x) 
%         = beta + [inv(x'x/n) ] var(x'e/n) [inv(x'x/n)]
% x'x/n         --p-->  M
% (n^-.5) x'e   --d-->  N(0, V)
%  =>var(x'e/n)=var(n^-.5 x'e)/n=V/n
%
%  var(b) = beta + M^-1 V/n M^-1
% V=X'(sigma^2 * I) X/n = sigma^2 M    (homoskedasticity)
%         = beta + M^-1 sigma^2/n
% (X' Md X)^-1 * (X' Md y)
%             Md * y = Md * X * b + resid 
% (X' Md X)^-1 * (X' Md * X * b + resid)
% b + (X' Md X)^ -1 resid

%XpMdXinv = inv(X'*Md*X );  %slow version
%XpMdXinv = inv(X'*MdX );


switch options.vce
  case 'ols'
    results.s2 = e' * e / (n - n_cat - k);
    results.bcov = inv(XprimMdX) * results.s2;
  
  case 'robust'
    % robust is equivalent to clustering on the category variable
    
    [V, n_c] = mdg_econometrics.getClusteredSMatrix(e, category_ids, X);
    results.bcov = (R \ (R' \ ((V / R) / R'))) * (n * ((n - 1) / (n - k - 1)) * (n_c / (n_c - 1))); 
        
    
  case 'cluster'    
    [V, n_c] = mdg_econometrics.getClusteredSMatrix(e, options.cluster_var, X);
    results.bcov = (R \ (R' \ ((V / R) / R'))) * (n * ((n - 1) / (n - k - 1)) * (n_c / (n_c - 1))); 
    % Notes:
    % The n is just to deal with that V is normalized by n
    %
    % The ((n - 1) / (n - k)) * (n_c / (n_c - 1))) is adjustment. Used by Stata, Hansen code etc...
    results.n_clusters = n_c;

  otherwise
    error(sprintf('Unknown VCE ''%s'' (should have already checked for this!?!?)',options.vce));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                 Fill out results struct                                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

results.bse = diag(results.bcov).^.5;
results.t = results.b ./ results.bse;

if(strcmp(options.vce, 'cluster') || strcmp(options.vce,'robust')) % both are basically clustering
    results.df_r = n_c - 1 ;           % residuals degree of freedom
else
    results.df_r = n - k - n_cat;
end
results.df_m = k;                      % model degrees of freedom

results.pvals = NaN(size(results.b));
results.pvals(results.t >= 0) = 2 * (1 - tcdf(results.t(results.t>=0), results.df_r));
results.pvals(results.t < 0)  = 2 * tcdf(results.t(results.t<0), results.df_r);
results.n_groups = n_cat;
results.n_obs    = n;

results.conf95  = [results.b + results.bse * tinv(.025, results.df_r), results.b + results.bse * tinv(.975, results.df_r)];

;
results.sigma_e = sqrt(e' * e / (n - n_cat - k));   % matches Stata


results.sigma_u = std(ubar);                         % matches Stata


results.F     = results.b' * inv(results.bcov) * results.b / k;
results.Fval = 1 - fcdf(results.F, results.df_m, results.df_r);





%u2 = u - mean(u);
%results.sigma_u2 = u2' * u2 ;

%asdf = corr(u, X * results.b)

%results.sigma_u = sqrt(u2' * u2 / (n - n_cat - k ));



%results.rss  = e'*e;                                   % residual sum of squares
%results.mss  = norm(yhat - mean(y), 2).^2;             % model sum of squares

% $$$ results.r2    = results.mss / (results.rss + results.mss);                           %n / var(y,1);

results.rss = e' * e;
results.r2_within  = corr(MdX * results.b, Mdy)^2;
results.r2_between = corr(Xbar * results.b, ybar)^2; 
results.r2_overall = corr(X * results.b, y)^2;
results.Tavg = n / n_cat;




% $$$ temp = ybari - mean(ybari);
% $$$ results.sigma_u = sqrt(temp' * temp / n);



%estimate sigma^2
