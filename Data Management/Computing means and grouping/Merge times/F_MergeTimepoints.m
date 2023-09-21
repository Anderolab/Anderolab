function [FreezTable] = F_MergeTimepoints(FreezTable, Groups, GroupNames)
%F_MERGETIMEPOINTS Summary of this function goes here
%   Detailed explanation goes here
for i = 1:length(Groups)
    GroupNames(i)
    Groups{i}
    FreezTable = F_CreateMeanVar(FreezTable, Groups{i}, GroupNames(i))
end
end

