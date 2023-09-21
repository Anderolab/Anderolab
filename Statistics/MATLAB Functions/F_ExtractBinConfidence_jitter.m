function [Mean, SD, Bonds] = F_ExtractBinConfidence_jitter(InputTable, ...
    Variable, BinWidth, Experiment, Iterations, Method)

% Determining size for randomising the indices
l = size(InputTable.(Variable), 2);
windowSize = 2700; % 90 seconds x 30 frames

% Convert End to frames
Experiment.FC.Task.End = Experiment.FC.Task.End * Experiment.FC.Task.FPS;

% Identify tone periods
toneIndices = find(contains(string(Experiment.FC.Task.Titles), 'Tone'));

% Addapting given the method parameter
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
    % Copy original data
    jitteredData = InputTable.(Variable);

    for idx = toneIndices
        % Calculate the start and end indices for the 90-second window around each tone
        startIdx = Experiment.FC.Task.Start(idx) - 900; % 30 seconds (or 900 frames) before the tone
        endIdx = Experiment.FC.Task.Start(idx) + 1800; % 60 seconds (or 1800 frames) including and after the tone

        % Extract the window
        windowData = jitteredData(:, startIdx:endIdx);
        
        % Jitterize within this window
        jitteredWindow = windowData(:, randperm(windowSize+1));
        % sizeLeft = size(jitteredData(:, startIdx:endIdx));
        % sizeRight = size(jitteredWindow);
        % 
        % disp('Left size:');
        % disp(sizeLeft);
        % disp('Right size:');
        % disp(sizeRight);

        % Replace the original window with the jitterized one
        jitteredData(:, startIdx:endIdx) = jitteredWindow;
    end
    InputTable.(Variable)=jitteredData;
    
    % Binning
    B(:, iter) = [mean(F_BinTraces(InputTable, Variable, ...
        BinWidth, Experiment), m); std(F_BinTraces(InputTable, Variable, ...
        BinWidth, Experiment), 0,m)];
    
    % Updating
    waitbar(iter/Iterations)
end
close(wb)
plot(B)

if strcmp(Method, 'Individualised')
    Mean = mean(mean(B(1:size(InputTable.(Variable), 1), :), m));
    SD = mean(mean(B(size(InputTable.(Variable), 1)+1:end, :), m));
else
    Mean = mean(B(1,:), m);
    %SD = std(B, [], m);
    SD = mean(B(2,:), m);
end
Bonds = [1, -1].*1.96.*SD + Mean;

end