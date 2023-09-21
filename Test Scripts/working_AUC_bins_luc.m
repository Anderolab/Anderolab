%% LOAD THE REQUIRED DATASETS
% 1 - Adding the project path to the searchpath
% Experiment properties
output_path = "C:\Users\Ander\OneDrive\Documents\MISCELL\TestEnv\Joaquín\TestFC";
% Adding functions to search path
addpath(genpath(output_path))

% 2 - Creating a new folder for the query
target_path = output_path + "\GlobalFluorescence Binning " + ...
    string(datetime(floor(now),'ConvertFrom','datenum'));
mkdir(target_path)
load("ExperimentData.mat")

%% 1 - Global fluorescence changes
% Extracting global flu
    % Setting the parameters for the specific query
    Data = Experiment.FC;
    Task = Experiment.FC.Task;
    
    % Extracting the palette and sexes
    GroupBy = string(Experiment.Project.Sexes);
    Palette = Experiment.Project.Palette;
    
    % Running the function
    Fluorescence = F_PopulationFluorescence(Data, GroupBy, "Sex", ...
        output_path);
    exportgraphics(gcf, strcat(target_path, "\Outliers.pdf"),  ...
        'ContentType','vector')

% Computing from mean fluorescence DeltaF/F0 - Shifting the traces
    
    % Parameters for the function
    RefEpoch = [];
    Method = "Pctile 50"; % "Mean", "Pctile n", "None"
    Scaling = true;
    
    % Actual shifting the traces
    Fluorescence = F_ShiftTrace(Fluorescence, 'Flu', [], ...
        Method, Task, Scaling);
%% 2 - Binning
% Actually binning
    % Parametters for the binning
    InputTable = Fluorescence;
    Variable = 'Flu';
    BinWidth = 5; %  In seconds

    % Actual binning
    [Binned, Repeated] = F_BinTraces(InputTable, Variable, BinWidth, ...
        Experiment);
        % Where binned are the mean global trace value for the specific bin
        % Frame-by-frame value of their corresponding bin
        % If we have 4-frame-long bins:
            % Binned = [1, 2, 3, 4]
            % Repeated = [1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4]

    % Saving
    Fluorescence.Binned = Binned;
    Fluorescence.RepBinned = Repeated;

% Extracting the baseline via randomisation
    Iterations = 1000;
    Method = "Global"; % "Individualised"
    [M, SD, Bonds] = F_ExtractBinConfidence(InputTable, Variable, BinWidth, ...
        Experiment, Iterations, Method);

%% 3 - Early data visualisation
% Visualising the traces
    F_ViewTraces(Fluorescence, 'RepBinned', 'Sex', Experiment, 30, 0, ...
        "Time (s)", "(\DeltaF/F_0)-F_{Hab}", output_path, ...
        false, [], false,[])
    yline(Bonds)
    exportgraphics(gcf, strcat(target_path, ...
        "\GlobalFluorescence_BIN.pdf"), 'ContentType','image') 

%% 4 - Grouping by tone
% Grouping
    % Inputs
    Task = Task;
    Episodes = ["Tone FC1", "Shock FC1";
        "Tone FC2", "Shock FC2"; ...
        "Tone FC3", "Shock FC3"
        "Tone FC4", "Shock FC4"; ...
        "Tone FC5", "Shock FC5"];
    Cathegory = ["Presentation (tone 1)", "Early (tones 2 & 3)", ...
        "Early (tones 2 & 3)", "Late (tones 4 & 5)", "Late (tones 4 & 5)"];
    Table = Fluorescence;
    Variable = 'RepBinned';
    Pre_Frames = 1000;
    Post_Frames = 2700;
    RF = 30;

    % Grouping function
    cat = F_MergeEpochs(Table, Variable, Episodes, Cathegory, Task, ...
        Pre_Frames, Post_Frames);

% Visualising
    SetZero = Pre_Frames;
    x_label = "Time (s)";
    y_label = "\DeltaF";
    close all

    F_ViewTraces(cat, 'RepBinned', ["Sex", "Cathegory"], Experiment, ...
        RF, SetZero, x_label, y_label, output_path, ...
        false, [], true, 29*30);

 % 5 - Computing confidence bonds given binned data (not shuffling)
% Computing mean +- 1.96SD
    mean(Fluorescence.Binned, "all")
    Z95 = 1.96*std(Fluorescence.Binned, [], "all");
    CI95 = mean(Fluorescence.Binned, "all") + [Z95, -Z95]

% Adding it to the figures
    ax = findall(gcf,'type','axes');
    add_interval = @(index) yline(ax(index), CI95(:).', "LineWidth", 2, ...
        'Color', 'r');
    for spt = 4:length(ax)
        add_interval(spt)
    end