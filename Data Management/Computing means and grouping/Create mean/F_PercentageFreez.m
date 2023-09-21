function [Freezing] = F_PercentageFreez(Behaviour, View)
%F_PERCENTAGEFREEZ Summary of this function goes here
%   Detailed explanation goes here
FreezPercentage = [];
Animals = string(fieldnames(Behaviour)).';

% Iterating per animal
for animal = Animals
    Epochs = string(fieldnames(Behaviour.(animal)));
    Epochs = Epochs(3:end);

    % Storage
    an_pctg = zeros(1, length(Epochs));
    
    % Iteration per epocj
    c = 1; % With a counter function
    for ep = Epochs.'
        % Computing %
        bh = Behaviour.(animal).(ep);
        an_pctg(c) = 100.*(sum(bh == 1)./length(bh));
        c = c+1;
    end

    % Storing
    FreezPercentage = [FreezPercentage; an_pctg];
    
end

Freezing = splitvars(table(FreezPercentage));
Freezing.Properties.VariableNames = Epochs;
Animal = Animals.';
Freezing = addvars(Freezing,Animal,'Before',Epochs(1));

if View == true
    imagesc(FreezPercentage)
    xticks(1:length(Epochs))
    xticklabels(Epochs)
    yticks(1:length(Animals))
    yticklabels(Animals)

    hold on
    yyaxis right
    plot(mean(FreezPercentage, 1), 'LineWidth', 3, "Color", 'w')
    colorbar
end
end

