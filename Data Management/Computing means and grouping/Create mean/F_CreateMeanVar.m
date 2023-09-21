function [FreezTable] = F_CreateMeanVar(FreezTable, Groups, Name)
%F_CREATEMEANVAR Summary of this function goes here
%   Detailed explanation goes here
FreezTable.(Name) = ...
    mean(table2array(FreezTable(:, ...
    ismember(string(FreezTable.Properties.VariableNames), Groups))), 2);
end

