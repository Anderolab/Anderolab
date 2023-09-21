function [Normalised] = F_95CI_Norm(Dataset)

Normalised = zeros(size(Dataset));
% 1 - Looping through the animals
for a_ix = 1:size(Dataset, 1)

    histogram(a_ix)
    hold on
    % 1 - Centering to the median
    Med_Shifted = Dataset(a_ix, :) - ...
        nanmedian(Dataset(a_ix, :));

    % 2 - Finding the confidence intervals
    CI = 1.96*nanstd(Med_Shifted);

    % 3 - Saving
    Normalised(a_ix, :) = Med_Shifted./CI;

end
hold off

