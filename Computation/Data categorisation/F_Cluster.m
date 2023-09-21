function [Clusters, K] = F_Cluster(DataMatrix)
%F_CLUSTER Summary of this function goes here
%   Detailed explanation goes here

%% Step 1 - Identifying the ideal number of clusters
% For speed, only selecting a number of datapoints
if size(DataMatrix, 1)>20000
    Ix = randperm(size(DataMatrix)); % Selection
    K = F_GetBestK(DataMatrix(Ix(1:30000), :), 200); % K identification
else
    K = F_GetBestK(DataMatrix, 200); % 
end



%% STEP 2 - Clustering
[Clusters] = kmeans(DataMatrix, K); % KMEANS
end

