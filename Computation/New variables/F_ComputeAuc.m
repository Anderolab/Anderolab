function [TableInput] = F_ComputeAuc(TableInput, Episodes, Cathegory, Task)
%F_COMPUTEAUC Summary of this function goes here
%   Detailed explanation goes here
AUC_s = zeros(size(TableInput, 1), length(unique(Cathegory)));
c_ = 1; % Counter
for cath_ = unique(Cathegory, "Stable")
    cath_
    eps_ = Episodes(Cathegory == cath_).'; % Identifying cath-specific eps.
    cath_AUCs = ...
        zeros(size(TableInput, 1), length(eps_)); % Cath save array
    c = 1; % Counter
    for ep_ = eps_
        s_ = Task.Start(string(Task.Titles) == ep_); % Start
        e_ = s_ + Task.Frames(string(Task.Titles) == ep_); % End
        cath_AUCs(:, c) = ...
            trapz(TableInput.Flu(:, s_:e_), 2); % Saving locally
        c = c+1; % Updating the counter
    end
    AUC_s(:, c_) = mean(cath_AUCs, 2);
    c_ = c_ +1;
end

% Saving
TableInput.AUC = AUC_s;

end

