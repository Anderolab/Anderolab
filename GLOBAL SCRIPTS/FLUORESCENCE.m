%% LOAD THE REQUIRED DATASETS
% 1 - Adding the project path to the searchpath
% Experiment properties
output_path = "C:\Users\Ander\OneDrive\Documents\MISCELL\TestEnv\Joaquín\TestFC";
% Adding functions to search path
addpath(genpath(output_path))

% 2 - Creating a new folder for the query
target_path = output_path + "\GlobalFluorescence final " + ...
    string(datetime(floor(now),'ConvertFrom','datenum'));
mkdir(target_path)
load("ExperimentData.mat")

%% 1 - Global fluorescence changes
% Setting the parameters for the specific query
Data = Experiment.FC;
Task = Experiment.FC.Task;
% Extracting the palette and sexes
GroupBy = string(Experiment.Project.Sexes);
Palette = Experiment.Project.Palette;
RF = 30;
PerformStats = false;
close all
% Running the function
Fluorescence = F_PopulationFluorescence(Data, GroupBy, "Sex", output_path);
exportgraphics(gcf, strcat(target_path, "\Outliers.pdf"),  ...
    'ContentType','vector')
% Shifting traces
RefEpoch = [];
Method = "Pctile 50"; % "Mean", "Pctile n", "None"
Scaling = true;

Fluorescence = F_ShiftTrace(Fluorescence, 'Flu', [], Method, Task, ...
      Scaling);

%%
close all
% Visualising the traces
F_ViewTraces(Fluorescence, 'Flu', 'Sex', Experiment, RF, 0, "Time (s)",...
    "(\DeltaF/F_0)-F_{Hab}", output_path, false, [], false,[])
exportgraphics(gcf, strcat(target_path, "\GlobalFluorescence.pdf"),  ...
    'ContentType','image')
close all


%% 3 - Grouping specific sessions
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
Variable = 'Flu';
Pre_Frames = 1000;
Post_Frames = 2700;
cat = F_MergeEpochs(Table, 'Flu', Episodes, Cathegory, Task, ...
    Pre_Frames, Post_Frames);

SetZero = Pre_Frames;
x_label = "Time (s)";
y_label = "\DeltaF";
close all
F_ViewTraces(cat, 'Flu', ["Sex", "Cathegory"], Experiment, RF, SetZero, ...
    x_label, y_label, output_path, false, [], true, 29*30);
exportgraphics(gcf, strcat(target_path, "\FC diffs.pdf"),  ...
    'ContentType','vector')

%% Computing area under the curve of each tone
Episodes = ["Tone FC1";
    "Tone FC2"; ...
    "Tone FC3";
    "Tone FC4"; ...
    "Tone FC5"];
Fluorescence_TEST = F_ComputeAuc(Fluorescence, Episodes, Cathegory, Task);
%%
close all
F_ViewTraces(Fluorescence_TEST, 'AUC', 'Sex', Experiment, 1, 0,...
    'Tone', '\DeltaF/F_0 (AUC)', char(target_path), true, '\StatsResults_AUC_Tone', ...
    false, []);
xticks([1:5])
box off
fig = gcf;
fig.Position = [861   287   230   443];
exportgraphics(gcf, strcat(target_path, "\FC AUC.pdf"),  ...
    'ContentType','vector')

Fluorescence

%% Computing area under the curve of each tone
close all
box off
Tones = zeros(size(Fluorescence, 1), 5);
for i = 1:5
    tone = strcat("Tone FC", num2str(i));
    frams = Task.Frames(string(Task.Titles) == tone);
    t_ = F_SelectEpoch(Fluorescence, 'Flu', tone, Task, 450, -frams);
    Tones(:, i) = trapz(t_.Flu, 2);
    
end
Fluorescence.Tones = Tones;
F_ViewTraces(Fluorescence, 'Tones', 'Sex', Experiment, 1, 0,...
    'Pre-tone ', '\DeltaF/F_0 (AUC)', char(target_path), true, ...
    '\StatsResults_AUC_Pretone', false, []);
xticks([1:5])
box off
fig = gcf;
fig.Position = [861   287   230   443];
exportgraphics(gcf, strcat(target_path, "\FC AUC_pretone.pdf"),  ...
    'ContentType','vector')
