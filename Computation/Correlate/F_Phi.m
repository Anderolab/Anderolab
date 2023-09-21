function [phi] = F_Phi(TP, TN, FP, FN)
%F_PHI Summary of this function goes here
%   Detailed explanation goes here
phi = (TP*TN - FN*FP)/sqrt((TP+ FP)*(TP + FN)*(FP + TN)*(FN + TN));
end

