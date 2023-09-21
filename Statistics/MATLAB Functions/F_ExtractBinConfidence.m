function [Mean, SD, Bonds] = F_ExtractBinConfidence(InputTable, ...
    Variable, BinWidth, Experiment, Iterations, Method)

% Determining size for randomising the indecies
l = size(InputTable.(Variable), 2);

% Addapting given the mehtod parameter
if Method == 'Individualised'
    B = zeros(2*size(InputTable.(Variable), 1), Iterations);
    m = 2;
else 
    B = zeros(2, Iterations);
    m = 'all';
end

% User update
wb = waitbar(0, "Extracting bin confidence interval");

% Iterating
for iter = 1:Iterations

    % Randomising
    InputTable.(Variable) = InputTable.(Variable)(:, randperm(l))
    
    % Binning
    B(:, iter) = [mean(F_BinTraces(InputTable, Variable, ...
        BinWidth, Experiment), m); std(F_BinTraces(InputTable, Variable, ...
        BinWidth, Experiment), 0,m)];
   

    % Updating
    waitbar(iter/Iterations)
end
close(wb)
plot(B)
% Computing outputs
if strcmp(Method, 'Individualised')
    Mean = mean(mean(B(1:size(InputTable.(Variable), 1), :), m));
    SD = mean(mean(B(size(InputTable.(Variable), 1)+1:end, :), m));
else
    Mean = mean(B(1,:), m);
    %SD = std(B, [], m);
    SD = mean(B(2,:), m);
end
Bonds = [1, -1].*1.96.*SD + Mean
end

