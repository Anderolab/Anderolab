function [Var, Distances, ClusterCentroids, ClustersSorted, NPerK] = ...
    F_VarPerK(Traces, Clusters)
%F_COMPUTEVARPERCLUSTER Summary of this function goes here
%   Detailed explanation goes here
ks = unique(Clusters);
ClusterCentroids = zeros(length(ks), size(Traces, 2));
Distances = [];
NPerK = [];
ClustersSorted = ks;
Var = zeros(length(ks), 1);
for k = 1:length(ks)
    ClusterCentroids(k, :) = mean(Traces(Clusters == ks(k), :), 1);
    distances_ = F_AssessDistances(Traces(Clusters == ks(k), :), ...
        ClusterCentroids(k, :));

    Var(k) = nanmean(distances_); %#ok<NANMEAN> 
    Distances = [Distances; distances_];
    NPerK = [NPerK, sum(Clusters == ks(k))];
        
end
end

