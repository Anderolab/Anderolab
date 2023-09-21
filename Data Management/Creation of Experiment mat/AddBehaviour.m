%% Setting
warning('off'); % Desabling warnings

% Behaviour properties
output_path = pwd;

% Expected animals
Animals = 1:14;
Sexes = {"Female", "Female", "Female", "Female", "Female", "Female", ...
    "Female", "Male", "Male", "Male", "Male", "Male", "Male", "Male"};
ColorDict = dictionary(["Male", "Female"], ...
    {[0, 75, 75]./255, [245, 173, 82]./255});

% Expected conditions
TrialTypes = ["FC", "FE1", "FE2"];


% Freez threshold
FreezThreshold = .8;

load('ExperimentData.mat')

%%
%% 1 - Finding all the relevant files
% Required function
    % F_ListFiles - Identifies all .mat and .xl files in the indicated
    % folder and its subfolders.

% Setting the file path
ms_path = "E:\Career\Raw Content\Work\Andero Lab - 2022-2023\" + ...
    "Computational\Calcium Imaging\Deconvoluted Raw Data\All Animals\Behav";

path_files = F_ListFiles(ms_path); % Identifying all the files in set path

% Creating the empty Behaviour output
Behaviour = struct([]);


%% 2 - Loading the files
% Report struct
Report = [];
Report.Missing = ["Missing file sessions:"]; %#ok<*NBRAK2>
Report.DeconvolutionError = ["Failed deconvolution sessions:"];
for trial_type = TrialTypes

    wb = waitbar(0, strcat("Saving behaviour for ", trial_type));
    
    % Gathering task-relevant information
    Task = [1, double(string(Experiment.(trial_type).Task.Lengths))];
    EpochTitles = string(Experiment.(trial_type).Task.Titles);
    start_ = cumsum(Task).';
    end_ = cumsum(Task(2:end)).';
    range = [start_(1:end-1), end_];

    % Looping through each animal
    for animal = Animals
        waitbar(animal/max(Animals))
        path_files
        % Identifying session-relevant path and saving it
        path_ = path_files(contains(path_files, ...
            F_NamingFunction_BH(animal, trial_type)));
        
        % Notifying and reporting missing files
        if isempty(path_) == 1
            prompt = strcat("    ", ...
                F_NamingFunction(animal, trial_type), " was not found.");
            fprintf('%s\n', prompt)
            Report.Missing = [Report.Missing; prompt];
            continue

        end

        % Saving the path
        Experiment.(trial_type).Behaviour.("M"+num2str(animal)).Path = ...
            path_;

        % Loading the behaviour
        BH = readtable(path_);
        BH = BH.Interval_3 > FreezThreshold;
        area(BH)

        % Saving the behaviour
        Experiment.(trial_type).Behaviour. ...
            ("M"+num2str(animal)).("AllFreez") = BH;

        % Cropping the behaviour
        for epoch = 1:size(range, 1)

            % Preventing code from stopping if the BH is not long enough
            try
                % Cropping and saving
                Experiment.(trial_type).Behaviour. ...
                    ("M"+num2str(animal)). ...
                    (regexprep(EpochTitles(epoch), ' ', '_')) = ...
                    BH(range(epoch, 1):range(epoch, 2));
            catch
                Experiment.(trial_type).Behaviour. ...
                    ("M"+num2str(animal)). ...
                    (regexprep(EpochTitles(epoch), ' ', '_')) = ...
                    "NO DATA";
            end
        end

        

    end
    close(wb)

end
%%
wb = waitbar(1, "Saving all data, please wait");
save('ExperimentBH.mat', "Experiment")
close(wb)