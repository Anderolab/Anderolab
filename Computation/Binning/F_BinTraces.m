function [Binned, Repeated] = F_BinTraces(InputTable, Variable, ...
    BinWidth, Experiment)

% Setting the parameters
RF = Experiment.Project.RF;
Width_Frams = RF*BinWidth;
N_Bins = floor(size(InputTable.(Variable), 2)/Width_Frams);

% Generating output
Binned = zeros(size(InputTable.(Variable), 1), N_Bins);

% Generating temporary start and end for identifying bin bounds
end_ = 0;
% Looping per bin
for bin = 1:N_Bins

    % Setting start and end
    start_ = end_ +1;
    end_ = start_ + Width_Frams - 1;

    % Computing mean flu for each neuron
    Binned(:, bin) = mean(InputTable.(Variable)(:, start_:end_), 2);
    
end

Repeated = repelem(Binned, 1, Width_Frams);
clear start_ end_ RF Width_Frams N_Bins
end

