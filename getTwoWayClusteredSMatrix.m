function [V, n_c1, n_c2] = getTwoWayClusteredSMatrix(resid, clusterids, X)
%
% This function is used behind the scenes by ivregress and others
%
% Requires: mdgtools.mg_getRowsWithKey
%
% [V, n_c] = getClusteredSMatrix(resid, clusterid, X)
%
%              resid: vector of residuals 
%          clusterids: cell array of two cluster ids or 
%                  X: right hand side matrix from regression
%
%                  V: the matrix we're estimating 
%                n_c: number of clusters

%%                                          START THE STUFF                                            %%
% c_unique = unique(c);                                     % METHOD 1: slow as shit and wasteful; use METHOD 2

%[n_rows, n_cols] = size(clusterids);

warning('two way cluster is not yet robustly tested');

if(iscell(clusterids))
    error('not yet supported (need to figure out the combined id stuff)');
    if(length(clusterids)==2)
        cluster_id1 = clusterids{1};
        cluster_id2 = clusterids{2};
    else
        if(size(clusterids,2)==2)
            cluster_id1 = clusterids(:,1);
            cluster_id2 = clusterids(:,2);
        else
            error('clusterids is cell but size of clusterids is funky');
        end
    end
else
    if(size(clusterids,2)==2)
        cluster_id1 = clusterids(:,1);
        cluster_id2 = clusterids(:,2);
        [~, ~, cluster_id12] = unique(clusterids,'rows');        
    else
        error('clusterids has bad size');
    end
end


[vce1, n_c1] = getClusteredSMatrix(resid, cluster_id1, X);
[vce2, n_c2] = getClusteredSMatrix(resid, cluster_id2, X);
[vce12, n_c12] = getClusteredSMatrix(resid, cluster_id12, X);

V = vce1 + vce2 - vce12;



% http://cameron.econ.ucdavis.edu/research/Cameron_Miller_Cluster_Robust_October152013.pdf