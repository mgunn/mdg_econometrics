function [results, options] = xtivregress(varargin)
% Estimate by instrumental variables
% Author: Matthew Gunn
% Note:   In general, behavior designed to match that of Stata
%
% [results, options] = ivregress(y, endo, exo, instruments, <option>, <option_val>, ...)
%
% Estimates y = [endo, exo] * b + epsilon%
%
% y:            N x 1
% endo:         N x L    endogenous variables we need to instrument
% exo :         N x K    A vector of 1s will be added to the end unless noconstant is set
% instruments:  N x L2   instruments for endogenous variables (L2 >= L)
%
% ids           N X 1    list of ids    
%
%    options:      option_values:
%
%      estimator:   This is the estimation technique used. Supported values are:
%                     '2sls' (equivalent to gmm with unadjusted)
%                     'gmm'  (default)
%                      
%        wmatrix:   If GMM estimation, these options are supported:
%                     'unadjusted'           Corresponds to 2sls estimator.
%                     'robust'  (default)
%                     'cluster'              If this is set, clustvar must also be set
%
%            vce:   By default, set to wmatrix type
%                     'unadjusted'
%                     'robust'
%                     'cluster' 
%
%    cluster_var:   Unique double or string to cluster on (required if clustering set)
%
%     noconstant:     1
%                     0  (default)
%
%  If more instruments
% 
% output:     results.b       estimate
%             results.bse     standard errors
%             results.bcov    covariance matrix

warning('MASSIVELY EXPERIMENTAL RIGHT NOW');    
    
%%                 Set Defaults                  %%
options = struct();
options.estimator = 'gmm';
options.noconstant = 0;

%%                 Set Options                   %%
if(nargin < 5)
    error('need at least 4 arguments: y, endo, exo, instruments, and ids');
end

y            = varargin{1};
endo         = varargin{2};
exo          = varargin{3};
instruments  = varargin{4};
category_ids = varargin{5};

if(isempty(exo))
    [y, endo, instruments, category_ids] = mdgtools.removenan(y, endo, instruments, category_ids);
else
    [y, endo, exo, instruments, category_ids] = mdgtools.removenan(y, endo, exo, instruments, category_ids);
end

i = 6;
while(i < nargin)
    switch varargin{i}
      case 'estimator'
        switch varargin{i+1}
          case '2sls'
            options.estimator = '2sls';
          case 'gmm'
            options.estimator = 'gmm';
          otherwise
            error(sprintf('unsupported estimator ''%s''', varargin{i+1}));
        end
      case 'vce'
        switch varargin{i+1}
          case 'unadjusted'
            options.vce = 'unadjusted';
          case 'robust'
            options.vce = 'robust';
          case 'cluster'
            options.vce = 'cluster';
          otherwise
            error(sprintf('Unknown vce ''%s''',varargin{i+1}));
        end
      case 'wmatrix'
        switch varargin{i+1}
          case 'unadjusted'
            options.wmatrix = 'unadjusted';
          case 'robust'
            options.wmatrix = 'robust';
          case 'cluster'
            options.wmatrix = 'cluster';
          otherwise
            error(sprintf('Unknown wmatrix ''%s''',varargin{i+1}));
        end
      case 'cluster_var'
        options.cluster_var = varargin{i+1};
      case 'noconstant'
        options.noconstant = 1;
      otherwise
        error(sprintf('Option ''%s'' is not understood', varargin{i}));        
    end %for switch
    i = i + 2;
end


%%                 Some defaults and checks               %%

% FIXME: check L2 >= L is OK... wasn't coded initially with that in mind

[N, cols] = size(y);

% $$$ if(~options.noconstant)
% $$$     exo = [exo, ones(N, 1)];
% $$$ end

if(cols > 1)
    error('only supports 1 dimensional y');
end;

[N2, L]  = size(endo);
[N3, K]  = size(exo);
[N4, L2] = size(instruments);
[N5, shoudlbe1] = size(category_ids);


%%                Set more defaults and checks             %%

switch options.estimator
  case 'gmm'
    if(~isfield(options,'wmatrix'))
        options.wmatrix = 'robust';     % Default is robust
    end
    if(~isfield(options,'vce'))
        options.vce = options.wmatrix;  % Default is wmatrix
    end

  case '2sls'
    if(~isfield(options,'vce'))
        options.vce = 'unadjusted';
    end
end

if(strcmp(options.vce,'cluster') || (isfield(options,'wmatrix') && strcmp(options.wmatrix,'cluster')))
    if(~isfield(options,'cluster_var'))
        error('options.cluster_var not set but clustering has been selected! Pass option ''cluster_var'' followed by cluster ids vector of cell array');
    end
end
    
    
%%                                     FINAL CHECK                                         %%
if((N ~= N2) || ((N2 ~=0) && (N2 ~= N)) || ((N3 ~= 0) && (N3 ~= N)) || ((N4 ~= 0) && (N4 ~= N)))
    error('Some input data has a number of rows inconsistent with the y vector');
end

%FIXME: check N = N2 = N3 = N4 and L = L2
if(L2 < L)
    error('Fewer instruments than endogenous regressors');
end

%%                      CONSTRUCT WITHIN AND BETWEEN DATA                       %%
% method 2: (idea: giant matrices are slower than using loops sensibly)
[uniqueid, map_seqid_to_catid] = mdgtools.mg_getRowsWithKey(category_ids);
n_cat = length(uniqueid);

%% within data  %%
MdXendo = zeros(N, L);                 % MdX: Whole X matrix with mean subtracted for each category_id
MdXexo  = zeros(N, K);                 % MdX: Whole X matrix with mean subtracted for each category_id
Mdy     = zeros(N, 1);                 % Mdy: Whole y vector with mean subtracted for each category_id

%% between data %%
ybar = zeros(n_cat, 1);            % Vector where each entry is the mean y for the category
Xbarendo = zeros(n_cat, L);            % Matrix where each entry is the mean X(i,:) for the category
Xbarexo  = zeros(n_cat, K);            % Matrix where each entry is the mean X(i,:) for the category

for i=1:length(uniqueid)
    list = map_seqid_to_catid{i};
    n_stu_i = length(list);

    Xtemp_endo = endo(list,:);
    Xtemp_exo  = exo(list,:);
    Xbar_endo_i = mean(Xtemp_endo);
    Xbar_exo_i = mean(Xtemp_exo);

    MdXendo(list,:) = Xtemp_endo - ones(n_stu_i,1) * Xbar_endo_i;
    MdXexo(list,:)  = Xtemp_exo  - ones(n_stu_i,1) * Xbar_exo_i;
    
    ytemp = y(list,:);
    ybar_i = mean(ytemp);
    Mdy(list) = ytemp - ybar_i;
    
    ybar(i)   = ybar_i;
    Xbarendo(i,:) = Xbar_endo_i;
    Xbarexo(i,:)  = Xbar_exo_i;
end



%%                              ESTIMATE results.b                                          %% 

X = [MdXendo, MdXexo];
Z = [MdXexo, instruments];

% $$$ X = [endo, exo];         % N x (L  + K)
% $$$ Z = [exo, instruments];  % N x (K + L2)

% These are the ingredients for all future linear algebra. Why do it more than once?
% Note: Could cholesky ZprimZ for repeated backsolve, but I don't think this is going to
% practically affect things much.

ZprimX = Z' * X;
Zprimy = Z' * y;
ZprimZ = Z' * Z;

if(L2 == L)  
    %System is exactly determined (unless colinear). Just solve linear system.
    results.b = linsolve(ZprimX, Zprimy);  % solve Z'y=Z'Xb  Note on dimension:  [(L2+K) x (L + K)] * x = [(L2 + K) x 1]
else 
    % In general, no unique solution if L2 > L. Use GMM weighting matrix    
    % Whether doing GMM or 2sls, we first need to do 2sls

% $$$     %% Method 1: 
% $$$     stage1_estimate = Z \ endo;
% $$$     endo_hat        = Z * stage1_estimate;   
% $$$     b_2sls          = [endo_hat, exo] \ y;
    

    %% Method 2: equivalent 
    % A = Z * inv(Z'*Z) * Z' * X
    % b_2sls = A \ y;
    b_2sls = (Z * (ZprimZ \ ZprimX)) \ y;
        
    if(strcmp(options.estimator,'2sls'))          % Note: corresponds to GMM with weighting matrix W = inv(Z'*Z);
        results.b = b_2sls;
    elseif(strcmp(options.estimator,'gmm'))
        % Step 1: First do two stage least squares
        u_2sls          = y - X * b_2sls;

        %% Basically, we calculate the S matrix and set W_inv = S
        if(strcmp(options.wmatrix,'unadjusted'))  
            W_inv = (u_2sls' * u_2sls / N) * ZprimZ / N;                    
        elseif(strcmp(options.wmatrix,'robust')) %do GMM two step procedure
            
            % Method 1 (matches Stata)            %calculate \frac{1}{n} \sum_i u_i^2 z_iz_i'
            temp = (u_2sls * ones(1, K + L2)) .* Z;
            W_inv = temp' * temp / N;
            
            % Method 2 (Doesn't match)
            % W_inv = cov((u_2sls * ones(1, K+L2)) .* Z,1)  % <----difference is that we're demeaning u_2sls
        elseif(strcmp(options.wmatrix,'cluster'))
            [W_inv, results.n_clusters] = mdg_econometrics.getClusteredSMatrix(u_2sls, options.cluster_var, Z);
        else
            error('unknown weighting matrix');
        end            


        
        %% GMM Deriviation for reference %%
        % Solve system: Z' * (y - X *b) = 0   <=> (Z'*X) * b = Z'y
        % Scenario 1: linear system
        % Scenario 2: L2 > L, i.e. overdetermined system
        % 
        % Then minimize (Z'*y - Z'*X*b)' * W * (Z'*y - Z'*X*b)
        %               (y'*Z - b'*X'*Z) * W * (Z'*y - Z'*X*b)
        %                y'ZWZ'y - b'X'ZWZ'y - y'ZWZ'Xb + b'X'ZWZ'Xb
        % b is only choice variable, W is symmetric therefore equivalent to minimize:
        %             b'X'ZWZ'Xb - y'ZWZ'Xb
        % Let A = X'*Z*W*Z'*X  Let c' = y'*Z*W*Z'*X        c = X'*Z * W * Z' * y
        % => minimize b'*A*b - c'*b                           %convex problem, FOC nec and suf
        % FOC: 2*A*b = c
                
        %see for more discussion : http://www.stata.com/manuals13/rivregress.pdf         
        % Now minimize the W norm (where W is pos def. weighting matrix)
        
        % This is what we want to do:
% $$$  %% ***************  METHOD 1 ***********
% $$$         W     = inv(W_inv);
% $$$         A     = ZprimX' * W * ZprimX;
% $$$         c     = ZprimX' * W * Zprimy;    % = (Zprimy' * W * ZprimX)'
% $$$         results.b = inv(A) * c;        
        
% $$$         %% METHOD 2:
% $$$         % This is a way without matrix inversion to do it. Solve linear systems rather than invert W_inv
        temp = W_inv \ ZprimX;          % Solve inv(W_inv) * ZprimX = W * ZprimX
        A    = ZprimX' * temp;
        c    = (Zprimy' * temp)';
        results.b = A \ c;        
        
% $$$         linsolve_opts = struct();
% $$$         linsolve_opts.POSDEF = true;
% $$$         linsolve_opts.SYM    = true;
% $$$         results.b = linsolve(A, c, linsolve_opts);
        
    else
        error('unknown estimator specified in options');
    end    
end


%%                                            just for convenient output                               %%
results.b_endo = results.b(1:L);
results.b_exo  = results.b(L+1:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                             Calculate covariance matrix of estimate                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


resid  = y - X * results.b;

%homoskedastic
%V = (resid'*resid/(N-(K+L))) * (Z'*Z/N);

%%FIXME: way too many matrix inversions etc.... this can be way better 

if(strcmp(options.estimator,'2sls'))
    W_inv = (resid'*resid / N) * (ZprimZ / N);
end

if(strcmp(options.vce,'unadjusted'))
    S = W_inv;
elseif(strcmp(options.vce,'robust'))
    temp = (resid * ones(1, K+L2)) .* Z;
    S = temp' * temp / N;
elseif(strcmp(options.vce,'cluster'))
    [S, results.n_clusters] = mdg_econometrics.getClusteredSMatrix(resid, options.cluster_var, Z);
else
    error('unknown vce type');
end
    
if(L == L2)               %%       exactly identified
    % What we want is:
    %         Minv = inv(M);
    %         results.bcov = Minv * V * Minv' / N;
    %
    % More precise way to numerically compute this is as follows (avoid matrix inversion)
    %    M = ZprimX / N;       
    %    results.bcov = (M \ (M \ S)')' / N;
    
    results.bcov = N * (ZprimX \ (S / ZprimX'));
else                      %%       overidentified situation 

    %% Algorithm 1:
% $$$     W     = inv(W_inv);
% $$$     M     = ZprimX / N;
% $$$     A     = M' * W * M;
% $$$     A_inv = inv(A);    
% $$$     results.bcov = A_inv * (M' * W * S * W * M) * A_inv / N;
% $$$     %results.bcov = inv(M'*W*M) * M' * W * S * W * M * inv(M'*W*M) / N;
% $$$     %             = 
    
    %% Algorithm 2:
    M      = ZprimX / N;
    WM     = W_inv \ M;
    MpWSWM = WM' * S * WM;    
    A      = M' * WM;
    results.bcov = (A \ MpWSWM) / A / N;
end


%%                                     STANDARD ERROR TIME                                 %%
results.bse = diag(results.bcov).^.5;
%results.r2 = NaN;

if(~options.noconstant)
    results.wald_chi2          = results.b(1:end-1)' * (results.bcov(1:end-1, 1:end-1) \ results.b(1:end-1));
end

results.z                     = results.b ./ results.bse;
results.pvals                 = NaN(size(results.b));
results.pvals(results.z >= 0) = 2 * (1 - normcdf(results.z(results.z>=0),0,1));
results.pvals(results.z < 0)  = 2 * normcdf(results.z(results.z<0),0,1);
results.conf95 = [results.b + results.bse * norminv(.025,0,1), results.b + results.bse * norminv(.975,0,1)];


%% Overidentification tests %%

%%FIXME: NEED TO OPEN POSSIBILIITY OF DIFFERENT W FOR ESTIMATION AND ERRORS

if(strcmp(options.estimator,'gmm') && (L2 > L)) % If we're doing GMM and have overidentifying restrictions

    %    W_zprimZ = inv(Z'*Z)

    %% Algorithm 1    
    
% $$$     W = inv(W_inv);
% $$$     temp = Zprimy - ZprimX * results.b;
% $$$     results.J = temp' * W * temp / N;

    %% Algorithm 2
    temp = Zprimy - ZprimX * results.b;
    results.J = temp' * (W_inv \ temp) / N; 
        
% $$$     %% Algorithm 3: wasteful and unnecessary
% $$$     Zprimresid     = (Z' * resid);
% $$$     results.J      = Zprimresid' * W * Zprimresid / N;

    results.J_chi2 = 1 - chi2cdf(results.J, L2 - L);
end

% $$$ results.t     = results.b ./ results.bse;
% $$$ results.df    = n - k - n_cat;
% $$$ results.pvals = NaN(size(results.b));
% $$$ results.pvals(results.t >= 0) = 2 * (1 - tcdf(results.t(results.t>=0), results.df));
% $$$ results.pvals(results.t < 0) = 2 * tcdf(results.t(results.t<0), results.df);



end
