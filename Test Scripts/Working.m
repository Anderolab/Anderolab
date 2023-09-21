Fluorescence
InputTable = Fluorescence;
Variable = 'Flu';
BinWidth = 5; %  In seconds

[Binned, Repeated] = F_BinTraces(InputTable, Variable, BinWidth, ...
    Experiment);

Iterations = 1000;
Method = "Global"; % "Individualised"
[M, SD, B] = F_ExtractBinConfidence(InputTable, Variable, BinWidth, ...
    Experiment, Iterations, Method);

M
SD
