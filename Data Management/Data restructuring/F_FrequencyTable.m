function [FrequencyTable] = F_FrequencyTable(A, B, Sort)
%F_FREQUENCYTABLE Summary of this function goes here
%   Detailed explanation goes here
A_ = reshape(unique(A, "stable"), 1, []);
B_ = reshape(unique(B, 'stable'), 1, []);

FrequencyTable = table();
FrequencyTable.G = B_.';

for a = A_
    freq = zeros(length(B_), 1);

    a_slice = B(A == a);
    c = 1;
    for b = B_
        freq(c) = sum(a_slice == b);
        c = c+1;
    end
    FrequencyTable.(string(a)) = (freq.*100)./sum(freq, 'all');
end

if Sort == true
    [~, ix] = sort(sum(table2array(FrequencyTable(:, 2:end)), 2), 'descend');
    FrequencyTable = FrequencyTable(ix, :);
end
end