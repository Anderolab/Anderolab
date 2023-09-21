function [Neuron] = F_SyntNeur(Times, FramPerPeak)

Neuron = [];

for t = 1:length(Times)
    ep_ = zeros(1, Times(t));
    ep_(1:floor(Times(t)/FramPerPeak(t))) = 1;
    ep_ = ep_(:, randperm(Times(t)));
    Neuron = [Neuron, ep_];
end
end

