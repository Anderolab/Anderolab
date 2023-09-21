%% Setting
warning('off'); % Desabling warnings

% Adding functions to search path
addpath(genpath("C:\Users\Ander\Desktop\Functions"))

% Experiment properties
output_path = uigetdir();

% Experiment name
ExpName = "TestLoad";

% Generating output file
output_path = strcat(output_path, '\', ExpName);
mkdir(output_path);
prompt = "Select where the animal configuration task mat file is found";
[file, path] = uigetfile('*.mat', prompt);
load(strcat(path, "\", file));


% Expected animals
Animals = length(AnimalConfig.Animals);
%Sexes = {"Female", "Female", "Female", "Female", "Female", "Female", ...
    %"Female", "Male", "Male", "Male", "Male", "Male", "Male", "Male"}; %#ok<CLARRSTR>

Groups=AnimalConfig.AnimalGroups;

ColorDict = dictionary(["Male", "Female"], ...
    {[0, 75, 75]./255, [245, 173, 82]./255});

% Expected conditions
%TrialTypes = ["FC", "FE1", "FE2"];
TrialTypes = ["IMO10", "IMO20"];

%RF
rf = 30;

%% 1 - Finding all the relevant files
% Required function
    % F_ListFiles - Identifies all .mat and .xl files in the indicated
    % folder and its subfolders.

% Setting the file path
ms_path = "C:\Users\Ander\Desktop\IMOsT2\Output 16 Aug 2023";
% addpath('C:\Users\Ander\OneDrive\Escritorio\Functions\File Management\File Search')
path_files = F_ListFiles(ms_path); % Identifying all the files in set path

% Creating the empty Experiment output
Experiment = struct([]);

%% 2 - Loading the files
% Report struct
Report = [];
Report.Missing = ["Missing file sessions:"];
Report.DeconvolutionError = ["Failed deconvolution sessions:"];
for trial_type = TrialTypes
%%

    % Populating the Output Struct
    Experiment(1).(trial_type) = [];

    % Looping through each animal
    for animal = Animals

        % Identifying session-relevant path and saving it
        path_ = path_files(contains(path_files, ...
            F_NamingFunction(animal, trial_type)));
        
        % Notifying and reporting missing files
        if isempty(path_) == 1
            prompt = strcat("    ", ...
                F_NamingFunction(animal, trial_type), " was not found.");
            fprintf('%s\n', prompt)
            Report.Missing = [Report.Missing; prompt];
            continue

        end

        % Saving the path
        Experiment.(trial_type).("M"+num2str(animal)).('Path') = path_;

        % Loading the ms file
        load(path_) % As MS

        % Verifying if the deconvolution was successfull
        % Notifying it, reporting it and removing it
        if sum(string(fieldnames(ms)) == "FiltTraces") == 0
            prompt = strcat("    ", ...
                F_NamingFunction(animal, trial_type), ...
                " failed during deconvolution.");
            fprintf('%s\n', prompt)
            Experiment.(trial_type) = rmfield(Experiment.(trial_type), ...
                "M"+num2str(animal));
            continue
        end

        
        raw_ = ms.RawTraces.';

        % Normalising with 5th percentile and saving         

        for n = 1:size(raw_, 1)
            raw_(n, :) = lowpass(raw_(n, :), 1, rf);
        end



        % Saving the raw trace
        Experiment.(trial_type).("M"+num2str(animal)).('Raw') = raw_;
        
        % Saving the filtered trace
        Experiment.(trial_type).("M"+num2str(animal)).('Filt') = ...
            ms.FiltTraces.';
        
    end
    % Asking the user if he wants to specifie the task file
    Variablein = input('If you want to add a task file press 1, if not press 2');
    if (Variablein == 1)
        % Asking the user to provide with task information
        prompt = "Select where the " + trial_type + ...
            " task mat file is found";
        [file, path] = uigetfile('*.mat', prompt); % TACKLE THIS
    
        % And adding it to the experiment dataset
        load(strcat(path, "\", file));
        Experiment.(trial_type).Task = Task;
        Experiment.(trial_type).Task.Location = ...
        strcat(path, "\", file);
    else 
        a='Task';
        Experiment.(trial_type).Task = a;
    end 
    
end
Report_ = ["STEP 1: LOADING"; ""; ""; Report.Missing; ""; ...
    Report.DeconvolutionError; ""];
clear Report
%% Saving
% Report
Experiment.Project = [];
Experiment.Project.Animals = Animals;
Experiment.Project.Groups = Groups;
Experiment.Project.Trials = TrialTypes;
Experiment.Project.RF = rf;
Experiment.Project.Sourcepath = ms_path;
Experiment.Project.Outputpath = output_path;
Experiment.Project.Palette = ColorDict;

writelines(Report_, strcat(output_path, "\Report.txt"));
save(strcat(output_path, "\ExperimentData.mat"), 'Experiment', '-v7.3');

