% 
% %% STEP 1 - LOADING THE REQUIRED DATASETS
output_path = "C:\Users\Ander\OneDrive\Documents\MISCELL\TestEnv\Joaquín\TestTunning\Pruebas validacion";
% 
% % Adding folder to search path
% addpath(genpath(output_path))
% 
% % 2 - Creating a new folder for the query
target_path = output_path + "\SingleNeuron " + ...
     string(datetime(floor(now),'ConvertFrom','datenum'));
mkdir(target_path)
% % 
% % 3 - Loading the required datasets
% load("ExperimentData.mat") % As Experiment, contains all traces
% load("SampleTask.mat") % As Task
%%
% STEP - SETTING THE INPUT FOR THE FUNCTION
% TunningFC = [];
Experiment.Project.Outputpath = output_path; % For testing purpuses
Neurons = Experiment.FE2;
Iterations = 10;
for i = 1:15
    Task = Experiment.FE2.Task;
    TunningEpochs = strcat("Tone FE", num2str(i));
    ReferenceEpochs = setdiff(string(Task.Titles), ...
        [TunningEpochs, "ITI FE15"]);
    [~, ~, ~,Output] = F_GetTunningProps_penal(Experiment, Task, ...
        TunningEpochs, ReferenceEpochs, Neurons, Iterations);
    TunningFC.(erase(TunningEpochs, ' ')) = Output;
end
save('TuningOutputFE2.mat', "TunningFC")

%%

% Extracting percentage of freez
Tones = fieldnames(TunningFE);
Tones = Tones(1:15)
Excited = [];
Inhibited = []; 
for tone = Tones.'

    % Finding column
    E_Field = TunningFE.(tone{:}).Properties.VariableNames(...
        contains(TunningFE.(tone{:}).Properties.VariableNames, ...
        ["Probability"]) + ...
        contains(TunningFE.(tone{:}).Properties.VariableNames, ...
        ["Excited"]) == 2)
    TunningFE.(tone{:})
    TunningFE.(tone{:}).(E_Field{:})
    Excited = [Excited, TunningFE.(tone{:}).(E_Field{:})]

    % Finding column
    I_Field = TunningFE.(tone{:}).Properties.VariableNames(...
        contains(TunningFE.(tone{:}).Properties.VariableNames, ...
        ["Probability"]) + ...
        contains(TunningFE.(tone{:}).Properties.VariableNames, ...
        ["Inhibited"]) == 2)
    TunningFE.(tone{:}).(I_Field{:})
    Inhibited = [Inhibited, TunningFE.(tone{:}).(I_Field{:})];
end

E_Grouped = [mean(Excited(:, [1:5]), 2), mean(Excited(:, [6:10]), 2),...
    mean(Excited(:, [11:15]), 2)]
I_Grouped = [Inhibited(:, [1:5]), mean(Inhibited(:, [6:10]), 2),...
    mean(Inhibited(:, [11:15]), 2)]

%% STATS 
% Tables Preparation Excited 
Excited = table();
Excited.Sex = TunningFE.(tone{:}).Sex;
Excited.Animal = TunningFE.(tone{:}).Animal;
Excited.Percentages = E_Grouped

writetable(Excited, "Excited.csv")

F_GlobalOutlierRemovalPlusSTATS('Excited.csv', 'Sex', 'Percentages', 'Animal', ...
    char(target_path), 'StatsResultsExcited')
F_ToneOutlierRemoval('Excited.csv', 'Sex', 'Percentages', 'Animal', ...
    char(target_path), 'StatsResultsExcited_ToneOutliers')
%%
% Tables Preparation Inhibited 
Inhibited = table();
Inhibited.Sex = TunningFE.(tone{:}).Sex;
Inhibited.Animal = TunningFE.(tone{:}).Animal;
Inhibited.Percentages = I_Grouped

writetable(Inhibited, "Inhibited.csv")

F_GlobalOutlierRemovalPlusSTATS('Inhibited.csv', 'Sex', 'Percentages', 'Animal', ...
    char(target_path), 'StatsResultsInhibited')
F_ToneOutlierRemoval('Inhibited.csv', 'Sex', 'Percentages', 'Animal', ...
    char(target_path), 'StatsResultsInhibited_ToneOutliers')
%%


E_means = [];
E_sds = [];
I_means = [];
I_sds = [];

% For visualisation
labels = ["Early FC", "Mid FC", "Late FC"]
yl = [];
% Extracting sex differences
Group = TunningFC.(tone{:}).Sex;
for g = unique(Group).'
    E_Grouped
    m_ = mean(E_Grouped(Group == g, :), 1)
    sd_ = std(E_Grouped(Group == g, :), [], 1)
    subplot(1, 2, 1)
    F_FillArea(m_, sd_./sqrt(sum(Group == g)), cell2mat(Experiment.Project.Palette(g)), 1:length(m_))
    hold on
    plot(m_, "LineWidth", 2, "Color", cell2mat(Experiment.Project.Palette(g)))
    hold on
    ylabel("% of responsive neurons")
    xticks(1:length(m_))
    xticklabels(labels)
    yl = [yl, ylim()];
    title("Excited")

    m_ = mean(I_Grouped(Group == g, :), 1)
    sd_ = std(I_Grouped(Group == g, :), [], 1)
    subplot(1, 2, 2)
    F_FillArea(m_, sd_./sqrt(sum(Group == g)), cell2mat(Experiment.Project.Palette(g)), 1:length(m_))
    hold on  
    plot(m_, "LineWidth", 2, "Color", cell2mat(Experiment.Project.Palette(g)))
    xticks(1:length(m_))
    xticklabels(labels)
    yl = [yl, ylim()];
    title("Inhibited")
end
ylim(get(gcf,'children'), [min(yl, [], 'all'), max(yl, [], 'all')])
%% 
% STEP - SETTING THE INPUT FOR THE MOVEMENT PROBABILITIES FUNCTION
Neurons = Experiment.IMO_L;
Iterations =  2;
[Output, FRECS] = F_GetMovProps(Experiment,...
     Neurons, Iterations)

%%
% STEP - SETTING THE INPUT FOR THE FUNCTION
TunningEpochs = ["Tone FC2", "Tone FC3"];
ReferenceEpochs = setdiff(string(Task.Titles), [TunningEpochs, "ITI FC5"]);
Neurons = Experiment.FC;


Output = F_GetTunningProps(Experiment, Task, TunningEpochs, ...
    ReferenceEpochs, Neurons, Iterations);
%%
% STEP - SETTING THE INPUT FOR THE FUNCTION
TunningEpochs = ["Tone FC1"];
ReferenceEpochs = setdiff(string(Task.Titles), [TunningEpochs, "ITI FC5"]);
Neurons = Experiment.FC;


Output = F_GetTunningProps(Experiment, Task, TunningEpochs, ...
    ReferenceEpochs, Neurons, Iterations);
 

%% STEP 1 - LOADING THE REQUIRED DATASETS
output_path = "C:\Users\Ander\OneDrive\Documents\MISCELL" + ...
    "\TestEnv\Joaquín\TestTunning";

% Adding folder to search path
addpath(genpath(output_path))

% 2 - Creating a new folder for the query
target_path = output_path + "\SingleNeuron " + ...
    string(datetime(floor(now),'ConvertFrom','datenum'));
mkdir(target_path)

% 3 - Loading the required datasets
load("ExperimentData.mat") % As Experiment, contains all traces

%% FOR CONDITIONING
% SETTING FUNCTION INPUTS
Experiment.Project.Outputpath="C:\Users\Ander\OneDrive\Documents\MISCELL" + ...
    "\TestEnv\Joaquín\TestTunning"
Neurons = Experiment.FC;
Iterations =  20;
Task = Experiment.FC.Task;

TunningEpochs = ["Tone FC4", "Tone FC5"];
ReferenceEpochs = setdiff(string(Task.Titles), [TunningEpochs, "ITI FC5"]);
%%
    % FOR THE PRESENTATION
    TunningEpochs = ["Pre-Tone FC1"];
    ReferenceEpochs = setdiff(string(Task.Titles), [TunningEpochs, "ITI FC5"]);
    %%
    Output = F_GetTunningProps(Experiment, Task, TunningEpochs, ...
        ReferenceEpochs, Neurons, Iterations);
%%
    % FOR EARLY CONDITIONING
    TunningEpochs = ["Pre-Tone FC2", "Pre-Tone FC3"];
    ReferenceEpochs = setdiff(string(Task.Titles), [TunningEpochs, "ITI FC5"]);
    
    
    Output = F_GetTunningProps(Experiment, Task, TunningEpochs, ...
        ReferenceEpochs, Neurons, Iterations);

    % FOR LATE CONDITIONING
    TunningEpochs = ["Pre-Tone FC4", "Pre-Tone FC5"];
    ReferenceEpochs = setdiff(string(Task.Titles), [TunningEpochs, "ITI FC5"]);
    
    
    Output = F_GetTunningProps(Experiment, Task, TunningEpochs, ...
        ReferenceEpochs, Neurons, Iterations);

%% FOR THE EXTINCTION 1
% SETTING FUNCTION INPUTS
Neurons = Experiment.FE1;
Iterations =  10;
Task = Experiment.FE1.Task;
    % FOR THE PRESENTATION
    TunningEpochs = ["Tone FE1"];
    ReferenceEpochs = setdiff(string(Task.Titles), [TunningEpochs, "ITI FE15"]);
    
    Output = F_GetTunningProps(Experiment, Task, TunningEpochs, ...
        ReferenceEpochs, Neurons, Iterations);

    % FOR EARLY CONDITIONING
    TunningEpochs = ["Tone FE2", "Tone FE3", "Tone FE4"];
    ReferenceEpochs = setdiff(string(Task.Titles), [TunningEpochs, "ITI FE5"]);
    
    
    Output = F_GetTunningProps(Experiment, Task, TunningEpochs, ...
        ReferenceEpochs, Neurons, Iterations);

    % FOR LATE CONDITIONING
    TunningEpochs = ["Tone FE13", "Tone FE14", "Tone FE15"];
    ReferenceEpochs = setdiff(string(Task.Titles), [TunningEpochs, "ITI FE5"]);
    
    
    Output = F_GetTunningProps(Experiment, Task, TunningEpochs, ...
        ReferenceEpochs, Neurons, Iterations);

%%
% STEP - SETTING THE INPUT FOR THE FUNCTION
TunningEpochs = ["Tone FC1"];
ReferenceEpochs = setdiff(string(Task.Titles), [TunningEpochs, "ITI FC5"]);

