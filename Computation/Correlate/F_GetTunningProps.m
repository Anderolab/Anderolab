
function [Output, Responses] = F_GetTunningProps(Experiment, Task, TunningEpochs, ...
    ReferenceEpochs, Dataset, Iterations)
%F_GETTUNNINGPROPS Summary of this function goes here
%   Detailed explanation goes here
%% STEP 1 - GENERATING THE STORAGE OUTPUTS
Expected = [];
Observed = [];
StandDevs = [];
STD_Distance = [];
Active = [];
Inhibited = [];
AllTraces = [];
Report = [];
lengths = [];
ResponseType = [];
AnimalPerNeuron = [];
SexPerNeuron = [];
OutputPath = Experiment.Project.Outputpath;
save_path = strcat(OutputPath, "\", join(TunningEpochs, " - "), ...
    " Tunning ", string(datetime(floor(now),'ConvertFrom','datenum'))); %#ok<TNOW1> 
prompt = strcat("   Results will be saved @ ", save_path);
Report = [Report; prompt; ""];
fprintf('% s\n', prompt)
% And the storage site
mkdir(save_path)


%% STEP 2 - IDENTIFYING OUTLIERS BEFORE PERFORMING THE TEST
% Identifying the animals
an_s = fieldnames(Dataset);
Animals = string(regexp(string(an_s(1:end-1)), '\d*', 'Match'));

% Finding all trial lengths;
c = 1;
for animal = Animals.'
    lengths(c) = size(Dataset.('M'+animal).("Raw"), 2);
    c = c+1;
end

% Identifying outliers
outliers = isoutlier(lengths, "mean");

% Reporting outliers to the user
for i = 1:sum(outliers)
    out_an = Animals(outliers);
    out_len = string(lengths(outliers));
    prompt = strcat("       Animal ", ...
        out_an(i), " was identified as an outlier with ",  ...
        out_len(i), ' frames.');
    fprintf('%s\n', prompt)
    Report = [Report; "   OUTLIER DETECTION:"; prompt];
end

% Visualising the outliers
boxplot(lengths)
xticks([])
ylabel("Frames (n)")

% Croppig all sessions and notifying the user
len = min(lengths);
prompt = strcat("       All sessions will be cropped to ", ...
    num2str(len), " frames.");
Report = [Report; prompt; ""];
fprintf('%s\n', prompt)

% Removing outliers
Animals = Animals(outliers == 0);

% Saving current figure
savename = strcat(save_path, "\Tunning - Outliers.pdf");
exportgraphics(gcf, savename, "ContentType", "vector")

%% STEP 3 - GENERATING MAIN USER'S AND STATISTICS OUTPUT TABLE
% Generating the output table
Output = table();
Output.Animal = double(Animals);
Sexes = string(Experiment.Project.Sexes).';
Output.Sex = Sexes(Output.Animal);

% And the variable names for the table storage
if length(TunningEpochs)>1
    Excited_ColName = strcat(join(TunningEpochs, " & "), " Excited IX");
    Inhibited_ColName = strcat(join(TunningEpochs, " & "), ...
        " Inhibited IX");
    ExcitedP_ColName = strcat("Probability of ", ...
        join(TunningEpochs, " & "), " Excited_1");
    InhibitedP_ColName = strcat("Probability of ", ...
        join(TunningEpochs, " & "), " Inhibited_2");
else 
    Excited_ColName = strcat(TunningEpochs, " Excited");
    Inhibited_ColName = strcat(TunningEpochs," Inhibited");
    ExcitedP_ColName = strcat("Probability of ", TunningEpochs, ...
        " Excited_1");
    InhibitedP_ColName = strcat("Probability of ", TunningEpochs, ...
        " Inhibited_2");
end
epoch={'Excibi,'}
% Generating the storage variable columns
Output.(Excited_ColName) = repelem({""}, length(Animals)).'; %#ok<STRSCALR> 
Output.(Inhibited_ColName) = repelem({""}, length(Animals)).'; %#ok<STRSCALR> 
Output.(ExcitedP_ColName) = zeros(length(Animals), 1);
Output.(InhibitedP_ColName) = zeros(length(Animals), 1);

% For figures
Ep_ = {ExcitedP_ColName, InhibitedP_ColName};
%% STEP 4 - GENERATING INDEXES FOR FRAMES OF INTEREST
% Identifying frames corresponding to the time of interest
TOI = []; % Time of interest
for Epoch = TunningEpochs
    ix = string(Task.Titles) == Epoch;
    TOI = [TOI, Task.Start(ix):(Task.Start(ix)+Task.Frames(ix)-30)];
end

% Identifying frames corresponding to the frames of reference
TOR = []; % Frames of reference
for Epoch = ReferenceEpochs
    ix = string(Task.Titles) == Epoch;
    TOR = [TOR, Task.Start(ix):(Task.Start(ix)+Task.Frames(ix))];
end

% Generating the binary for the comparison
BinarisedTask = repelem("0", len);
BinarisedTask(TOI) = "1";
prompt = "   Threshold will be of 1.96 SDs from mean (95% CI)";
Report = [Report; ""; prompt];
fprintf('%s\n', prompt)
%% STEP 5 - PERFORMING THE TEST
% Looping through animals
c = 1; % Counter function
for an = Animals.'
    wb_ = waitbar(0, ...
        strcat("Identifying responsive neurons in animal ", num2str(an)));
    prompt = strcat("   Processing animal ", num2str(an));
    fprintf('%s\n', prompt)
    Report = [Report; prompt];
    % Gathering the animal specific data
    Maxims = double(islocalmax(Dataset.(strcat('M', ...
        num2str(an))).Filt(:, 1:len), 2));
    Peaks = string(Maxims);
    Intersect = BinarisedTask+Peaks;
    % Identifying the parameters for each neuron
    oo = sum(Intersect == "00", 2);
    lo = sum(Intersect == "10", 2);
    ol = sum(Intersect == "01", 2);
    ll = sum(Intersect == "11", 2);
    
    Obv = F_ComputePhi(oo, ol, lo, ll);
    % Computing Phi
    Observed = [Observed; Obv];

    % Iterating to attain the expected
    Expected_Scores = gpuArray(zeros(size(Maxims, 1), Iterations));

    for iter = 1:Iterations

        waitbar(iter/Iterations);
        Rand_Peaks = Peaks(:, randperm(len));
        Intersect_Rand = BinarisedTask+Rand_Peaks;
        oo = sum(Intersect_Rand == "00", 2);
        lo = sum(Intersect_Rand == "10", 2);
        ol = sum(Intersect_Rand == "01", 2);
        ll = sum(Intersect_Rand == "11", 2);
        Expected_Scores(:, iter) = F_ComputePhi(oo, ol, lo, ll);

    end
    Exp = mean(Expected_Scores, 2);
    Expected = [Expected; Exp];
    SDs = std(Expected_Scores, [], 2);
    StandDevs = [StandDevs; SDs];
    Distances = (Obv-Exp)./SDs;
    STD_Distance = [STD_Distance; Distances];

    % Saving the activated and inhibited neurons
    Active = [Active; ...
        Dataset.(strcat('M', num2str(an))).Filt(Distances > 1.96, 1:len)];
    Inhibited = [Inhibited; ...
        Dataset.(strcat('M', num2str(an))).Filt(Distances < -1.96, 1:len)];
    % All traces
    AllTraces = [AllTraces; ...
        Dataset.(strcat('M', num2str(an))).Filt(:, 1:len)];

    % Saving
    Excit = find(Distances > 1.96);
    Inhibit = find(Distances < -1.96);
    Output.(Excited_ColName)(c) = {Excit};
    Output.(Inhibited_ColName)(c) = {Inhibit};
    Output.(ExcitedP_ColName)(c) = 100*length(Excit)/size(Maxims, 1);
    Output.(InhibitedP_ColName)(c) = 100*length(Inhibit)/size(Maxims, 1);
    
        % For the frequency test
        Tunning_ = repelem("Unresponsive", length(Distances));
        Tunning_(Excit) = "Excited";
        Tunning_(Inhibit) = "Inhibited";
        ResponseType = [ResponseType, Tunning_];

        % Saving the animal
        AnimalPerNeuron = [AnimalPerNeuron, ...
            repelem(an, length(Distances))];
        SexPerNeuron = [SexPerNeuron, ...
            repelem(Experiment.Project.Sexes{double(an)}, ...
            length(Distances))];

    prompt = strcat("       ",  num2str(length(Excit)), ...
        " stimulus-excited neurons have been identified for animal ", ...
        num2str(an));
    fprintf('%s\n', prompt)
    Report = [Report; prompt];
    prompt = strcat("       ",  num2str(length(Inhibit)), ...
        " stimulus-inhibited neurons have been identified for animal ", ...
        num2str(an));
    fprintf('%s\n', prompt)
    Report = [Report; prompt];    
    c = c +1;
    close(wb_);

end

%% STEP 6 - GENERATING THE VISUALISATIONS
% First figure - Methods
    [~, sort_ix] = sort(Observed);
    F_FillArea(Expected(sort_ix).', (StandDevs(sort_ix).*1.96).', ...
        'k', 1:length(Expected(sort_ix)))
    hold on
    plot(Expected(sort_ix), "Color", 'k')
    hold on
    plot(Observed(sort_ix), "Color", 'r', "LineWidth", 2)
    O = Observed(sort_ix);
    STD_Sorted = STD_Distance(sort_ix);
    sig_ix = find(abs(STD_Sorted) > 1.96);
    scatter(sig_ix, O(sig_ix), 20, 'K', 'filled')
    hold off
    legend(["95% CI", "Expected", "Observed", "Significant"], ...
        "Location","northwest");
    xlim([1, length(Observed)])
    xlabel("Neuron", "FontSize", 12)
    ylabel("\phi Coefficient", "FontSize", 12)
    set(gcf,'Position',[400 100 300 400])
    box off
    hold off
    savename = strcat(save_path, "\Tunning - Neuron identification.pdf");
    exportgraphics(gcf, savename, "ContentType","vector")

    
% Second figure - Sample
    % For the active neurons - Calculating the error
    close all
    sem_ = std(Active, [], 1)./sqrt(size(Active, 1));
    
    % Getting the axes
    F_FillArea(mean(Active, 1), sem_, [212, 100, 66]./255, ...
        1:length(Active))
    yl = ylim();
    
    % Labelling areas of interest
    area(double(BinarisedTask), 'EdgeColor','none', "FaceAlpha", .3, ...
        'FaceColor', 'k')
    hold on

    % Visualising the errors for the active neurons
    F_FillArea(mean(Active, 1), sem_, [212, 100, 66]./255, ...
        1:length(Active))
    % Viewing active neurons
    plot(mean(Active, 1), "Color", [212, 100, 66]./255)
    hold on
    F_FillArea(mean(Inhibited, 1), sem_, [76, 148, 199]./255, ...
        1:length(Inhibited))
    plot(mean(Inhibited, 1), "Color", [76, 148, 199]./255)
    
    ylim([0, max(yl)])
    xlim([1, length(Active)])
    set(gcf, 'Position', [400, 100, 1200, 500])
    box off
    xlabel("Time (Frames)")
    ylabel("Mean Fluorescence")
    legend([join(TunningEpochs, ' & '), "", "Excited", "", "Inhibited"])
    hold off
    savename = strcat(save_path, "\Tunning - Group mean.pdf");
    exportgraphics(gcf, savename, "Resolution", 1200)

% Third figure - Individual-neuron level visualisation
    close all
    % Selecting the top five neurons;
    Sorted = sort(Observed);
    TopInactive = AllTraces(Observed < Sorted(6), :);

    % Getting the axes
    area(double(BinarisedTask).*13, 'EdgeColor','none', ...
        "FaceAlpha", .3, 'FaceColor', 'k')
    hold on
    for i = 1:5
        plot(TopInactive(i, :) + i, "Color", [76, 148, 199]./255)
        hold on
    end
    Sorted = flip(Sorted);
    TopActive = AllTraces(Observed > Sorted(6), :);
    for i = 7:11

        plot(TopActive(i-6, :)+i, "Color", [212, 100, 66]./255)
        hold on
    end
    hold off
    ylim([0, 13])
    yticks([3, 9])
    yticklabels({"Inhibited", "Excited"}); %#ok<CLARRSTR> 
    ytickangle(90)
    xlabel("Time (Frames)")
    ylabel("Filtered fluorescence")
    xlim([1, size(AllTraces, 2)])
    box off
    Fig_ = gcf;
    Fig_.Position = [400, 100, 600, 500];
    savename = strcat(save_path, "\Tunning - Sample Neurons.pdf");
    exportgraphics(gcf, savename, "ContentType", "vector")

% Fourth figure - Ratios
    % Identifying the percentage of neuron tunning
    close all
    Properties = [sum(Distances > 1.96, 'All'), ...
        sum(Distances < -1.96, 'All')];
    Properties(end+1) = length(Distances) - sum(Properties);
    pie(Properties)
    ax = gca();
    ax.Colormap = [212, 100, 66; 76, 148, 199; 200, 200, 200]./255; 
    legend(["Excited", "Inhibited", "Irresponsive"], "Location", ...
        "northeast")
    F_ = gcf;
    F_.Position = [400, 100, 390, 320];
    savename = strcat(save_path, "\Tunning - Ratios.pdf");
    exportgraphics(gcf, savename, "ContentType", "vector")
 
% Fifth figure - Distinction between conditions
    % Generating the first subplot
    close all
    yls = [];
    subplot(1, 2, 1)

    boxplot(Output.(ExcitedP_ColName), Output.Sex, ...
        "Colors", 'k');
    colours = ...
        flip(cell2mat(Experiment.Project.Palette(unique(Output.Sex, ...
        'stable'))), 1);
    f = findobj(gca, 'Tag', 'Box');
    for i = 1:length(f)
        patch(get(f(i),'XData'),get(f(i),'YData'),colours(i,:), ...
            'FaceAlpha',.6);
    end
    h=findobj('LineStyle','--'); set(h, 'LineStyle','-');
    box off
    hold on
    ylabel("Tunned neurons (%)")
    xlabel("Excited")
    yls = [yls; ylim()];

    subplot(1, 2, 2)
    boxplot(Output.(InhibitedP_ColName), Output.Sex, ...
        "Colors", 'k');
    f = findobj(gca, 'Tag', 'Box');
    for i = 1:length(f)
        patch(get(f(i),'XData'),get(f(i),'YData'),colours(i,:), ...
            'FaceAlpha',.6);
    end
    h=findobj('LineStyle','--'); set(h, 'LineStyle','-');
    xlabel("Inhibited")
    h = gca;
    h.YAxis.Visible = 'off';
    box off
    yls = [yls; ylim()];
    ylim([min(yls(:, 1)), max(yls(:, 2))])

    
    savename = strcat(save_path, "\Tunning - Group Differences.pdf");
    exportgraphics(gcf, savename, "ContentType", "vector")
    close all

%% STEP 7 - CLOSING
% Saving the report
savename = strcat(save_path, "\Report.txt");
writelines(Report,savename);
% And the data
savename = strcat(save_path, "\TunningOutput.mat");
Tunning = Output;
save(savename, "Tunning");

% Generating the frequency table
Responses = table();
size(SexPerNeuron)
size(ResponseType)
size(AnimalPerNeuron)
Responses.Sex = SexPerNeuron.';
Responses.Animal = AnimalPerNeuron.';
Responses.ResponseType = ResponseType.';

% And saving it
savename = strcat(save_path, "\SingleNeuronResponses.mat");
save(savename, "Responses");
savename_chi = strcat(save_path, "\SingleNeuronResponses.csv");
writetable(Responses, savename_chi);

%% STEP 8 - STATS
% Save Table of Data
savename = strcat(save_path, "\TunningOutput.csv");
StatsTable = table();
StatsTable.Sex = Output.Sex;
StatsTable.Animal = Output.Animal;
StatsTable.(ExcitedP_ColName) = Output.(ExcitedP_ColName);
StatsTable.(InhibitedP_ColName) = Output.(InhibitedP_ColName);
writetable(StatsTable, savename);

% Make Table of Variables 
GroupBy = "Sex";
Variable1 = ExcitedP_ColName;
Variable2 = InhibitedP_ColName;
Id = "Animal";
DatasetPath = string(savename);
NewFolderName = "StatResults";
NewFolderPath = strcat(save_path, "\", NewFolderName);
DatasetPath_ChiTest = string(savename_chi)

F_MakeCSVForStat(GroupBy, Variable1, Variable2, Id, DatasetPath, NewFolderPath, DatasetPath_ChiTest);

% Run R Script of Analysis
ScriptPath = strcat('', pwd, '\BoxTestAndANOVA.R')
F_RunRScript(ScriptPath)


end


% function [Output, Responses] = F_GetTunningProps(Experiment, Task, TunningEpochs, ...
%     ReferenceEpochs, Dataset, Iterations)
% %F_GETTUNNINGPROPS Summary of this function goes here
% %   Detailed explanation goes here
% %% STEP 1 - GENERATING THE STORAGE OUTPUTS
% Expected = []; % From shuffling
% Observed = []; % From original
% StandDevs = []; % From shuffling
% STD_Distance = []; % Distance in STDs
% Active = []; % Active neuron traces (original)
% Inhibited = []; % Inhibited neuron traces (original)
% AllTraces = []; % All traces (original)
% Report = []; % Reporting the user of the progress
% lengths = []; % Lenths for all sessions
% OutputPath = Experiment.Project.Outputpath; % Gathering the output path
%                         % Saved in the experiment struct
% 
% % For the output table
% ResponseType = [];
% AnimalPerNeuron = [];
% SexPerNeuron = [];
% save_path = strcat(OutputPath, "\", join(TunningEpochs, " - "), ...
%     " Tunning ", string(datetime(floor(now),'ConvertFrom','datenum'))); %#ok<TNOW1> 
% prompt = strcat("   Results will be saved @ ", save_path);
% Report = [Report; prompt; ""];
% fprintf('%s\n', prompt)
% % And the storage site
% mkdir(save_path)
% 
% 
% %% STEP 2 - IDENTIFYING OUTLIERS BEFORE PERFORMING THE TEST
% % Identifying the animals
% an_s = fieldnames(Dataset);
% Animals = string(regexp(string(an_s(1:end-1)), '\d*', 'Match'));
% 
% % Finding all trial lengths;
% c = 1;
% for animal = Animals.'
%     lengths(c) = size(Dataset.('M'+animal).("Raw"), 2);
%     c = c+1;
% end
% 
% % Identifying outliers
% outliers = isoutlier(lengths, "mean");
% 
% % Reporting outliers to the user
% for i = 1:sum(outliers)
%     out_an = Animals(outliers);
%     out_len = string(lengths(outliers));
%     prompt = strcat("       Animal ", ...
%         out_an(i), " was identified as an outlier with ",  ...
%         out_len(i), ' frames.');
%     fprintf('%s\n', prompt)
%     Report = [Report; "   OUTLIER DETECTION:"; prompt];
% end
% 
% % Visualising the outliers
% boxplot(lengths)
% xticks([])
% ylabel("Frames (n)")
% 
% % Croppig all sessions and notifying the user
% len = min(lengths);
% prompt = strcat("       All sessions will be cropped to ", ...
%     num2str(len), " frames.");
% Report = [Report; prompt; ""];
% fprintf('%s\n', prompt)
% 
% % Removing outliers
% Animals = Animals(outliers == 0);
% 
% % Saving current figure
% savename = strcat(save_path, "\Tunning - Outliers.pdf");
% exportgraphics(gcf, savename, "ContentType", "vector")
% 
% %% STEP 3 - GENERATING MAIN USER'S AND STATISTICS OUTPUT TABLE
% % Generating the output table
% Output = table();
% Output.Animal = double(Animals);
% Sexes = string(Experiment.Project.Sexes).';
% Output.Sex = Sexes(Output.Animal);
% 
% % And the variable names for the table storage
% if length(TunningEpochs)>1
%     Excited_ColName = strcat(join(TunningEpochs, " & "), " Excited IX");
%     Inhibited_ColName = strcat(join(TunningEpochs, " & "), ...
%         " Inhibited IX");
%     ExcitedP_ColName = strcat("Probability of ", ...
%         join(TunningEpochs, " & "), " Excited_1");
%     InhibitedP_ColName = strcat("Probability of ", ...
%         join(TunningEpochs, " & "), " Inhibited_2");
% else 
%     Excited_ColName = strcat(TunningEpochs, " Excited");
%     Inhibited_ColName = strcat(TunningEpochs," Inhibited");
%     ExcitedP_ColName = strcat("Probability of ", TunningEpochs, ...
%         " Excited_1");
%     InhibitedP_ColName = strcat("Probability of ", TunningEpochs, ...
%         " Inhibited_2");
% end
% 
% % Generating the storage variable columns
% Output.(Excited_ColName) = repelem({""}, length(Animals)).'; %#ok<STRSCALR> 
% Output.(Inhibited_ColName) = repelem({""}, length(Animals)).'; %#ok<STRSCALR> 
% Output.(ExcitedP_ColName) = zeros(length(Animals), 1);
% Output.(InhibitedP_ColName) = zeros(length(Animals), 1);
% 
% % For figures
% Ep_ = {ExcitedP_ColName, InhibitedP_ColName};
% %% STEP 4 - GENERATING INDEXES FOR FRAMES OF INTEREST
% % Identifying frames corresponding to the time of interest
% TOI = []; % Time of interest
% for Epoch = TunningEpochs
%     ix = string(Task.Titles) == Epoch;
%     TOI = [TOI, Task.Start(ix):(Task.Start(ix)+Task.Frames(ix)-30)];
% end
% 
% % Identifying frames corresponding to the frames of reference
% TOR = []; % Frames of reference
% for Epoch = ReferenceEpochs
%     ix = string(Task.Titles) == Epoch;
%     TOR = [TOR, Task.Start(ix):(Task.Start(ix)+Task.Frames(ix))];
% end
% 
% % Generating the binary for the comparison
% BinarisedTask = repelem("0", len);
% BinarisedTask(TOI) = "1";
% prompt = "   Threshold will be of 1.96 SDs from mean (95% CI)";
% Report = [Report; ""; prompt];
% fprintf('%s\n', prompt)
% %% STEP 5 - PERFORMING THE TEST
% % Looping through animals
% c = 1; % Counter function
% for an = Animals.'
%     wb_ = waitbar(0, ...
%         strcat("Identifying responsive neurons in animal ", num2str(an)));
%     prompt = strcat("   Processing animal ", num2str(an));
%     fprintf('%s\n', prompt)
%     Report = [Report; prompt];
%     % Gathering the animal specific data
%     Maxims = double(islocalmax(Dataset.(strcat('M', ...
%         num2str(an))).Filt(:, 1:len), 2));
%     Peaks = string(Maxims);
%     Intersect = BinarisedTask+Peaks;
%     % Identifying the parameters for each neuron
%     oo = sum(Intersect == "00", 2);
%     lo = sum(Intersect == "10", 2);
%     ol = sum(Intersect == "01", 2);
%     ll = sum(Intersect == "11", 2);
% 
%     Obv = F_ComputePhi(oo, ol, lo, ll);
%     % Computing Phi
%     Observed = [Observed; Obv];
% 
%     % Iterating to attain the expected
%     Expected_Scores = gpuArray(zeros(size(Maxims, 1), Iterations));
% 
%     for iter = 1:Iterations
% 
%         waitbar(iter/Iterations);
%         Rand_Peaks = Peaks(:, randperm(len));
%         Intersect_Rand = BinarisedTask+Rand_Peaks;
%         oo = sum(Intersect_Rand == "00", 2);
%         lo = sum(Intersect_Rand == "10", 2);
%         ol = sum(Intersect_Rand == "01", 2);
%         ll = sum(Intersect_Rand == "11", 2);
%         Expected_Scores(:, iter) = F_ComputePhi(oo, ol, lo, ll);
% 
%     end
%     Exp = mean(Expected_Scores, 2);
%     Expected = [Expected; Exp];
%     SDs = std(Expected_Scores, [], 2);
%     StandDevs = [StandDevs; SDs];
%     Distances = (Obv-Exp)./SDs;
%     STD_Distance = [STD_Distance; Distances];
% 
%     % Saving the activated and inhibited neurons
%     Active = [Active; ...
%         Dataset.(strcat('M', num2str(an))).Filt(Distances > 1.96, 1:len)];
%     Inhibited = [Inhibited; ...
%         Dataset.(strcat('M', num2str(an))).Filt(Distances < -1.96, 1:len)];
%     % All traces
%     AllTraces = [AllTraces; ...
%         Dataset.(strcat('M', num2str(an))).Filt(:, 1:len)];
% 
%     % Saving
%     Excit = find(Distances > 1.96);
%     Inhibit = find(Distances < -1.96);
%     Output.(Excited_ColName)(c) = {Excit};
%     Output.(Inhibited_ColName)(c) = {Inhibit};
%     Output.(ExcitedP_ColName)(c) = 100*length(Excit)/size(Maxims, 1);
%     Output.(InhibitedP_ColName)(c) = 100*length(Inhibit)/size(Maxims, 1);
% 
%         % For the frequency test
%         Tunning_ = repelem("Unresponsive", length(Distances));
%         Tunning_(Excit) = "Excited";
%         Tunning_(Inhibit) = "Inhibited";
%         ResponseType = [ResponseType, Tunning_];
% 
%         % Saving the animal
%         AnimalPerNeuron = [AnimalPerNeuron, ...
%             repelem(an, length(Distances))];
%         SexPerNeuron = [SexPerNeuron, ...
%             repelem(Experiment.Project.Sexes{double(an)}, ...
%             length(Distances))];
% 
%     prompt = strcat("       ",  num2str(length(Excit)), ...
%         " stimulus-excited neurons have been identified for animal ", ...
%         num2str(an));
%     fprintf('%s\n', prompt)
%     Report = [Report; prompt];
%     prompt = strcat("       ",  num2str(length(Inhibit)), ...
%         " stimulus-inhibited neurons have been identified for animal ", ...
%         num2str(an));
%     fprintf('%s\n', prompt)
%     Report = [Report; prompt];    
%     c = c +1;
%     close(wb_);
% 
% end
% 
% %% STEP 6 - GENERATING THE VISUALISATIONS
% % First figure - Methods
%     [~, sort_ix] = sort(Observed);
%     F_FillArea(Expected(sort_ix).', (StandDevs(sort_ix).*1.96).', ...
%         'k', 1:length(Expected(sort_ix)))
%     hold on
%     plot(Expected(sort_ix), "Color", 'k')
%     hold on
%     plot(Observed(sort_ix), "Color", 'r', "LineWidth", 2)
%     O = Observed(sort_ix);
%     STD_Sorted = STD_Distance(sort_ix);
%     sig_ix = find(abs(STD_Sorted) > 1.96);
%     scatter(sig_ix, O(sig_ix), 20, 'K', 'filled')
%     hold off
%     legend(["95% CI", "Expected", "Observed", "Significant"], ...
%         "Location","northwest");
%     xlim([1, length(Observed)])
%     xlabel("Neuron", "FontSize", 12)
%     ylabel("\phi Coefficient", "FontSize", 12)
%     set(gcf,'Position',[400 100 300 400])
%     box off
%     hold off
%     savename = strcat(save_path, "\Tunning - Neuron identification.pdf");
%     exportgraphics(gcf, savename, "ContentType","vector")
% 
% 
% % Second figure - Sample
%     % For the active neurons - Calculating the error
%     close all
%     sem_ = std(Active, [], 1)./sqrt(size(Active, 1));
% 
%     % Getting the axes
%     F_FillArea(mean(Active, 1), sem_, [212, 100, 66]./255, ...
%         1:length(Active))
%     yl = ylim();
% 
%     % Labelling areas of interest
%     area(double(BinarisedTask), 'EdgeColor','none', "FaceAlpha", .3, ...
%         'FaceColor', 'k')
%     hold on
% 
%     % Visualising the errors for the active neurons
%     F_FillArea(mean(Active, 1), sem_, [212, 100, 66]./255, ...
%         1:length(Active))
%     % Viewing active neurons
%     plot(mean(Active, 1), "Color", [212, 100, 66]./255)
%     hold on
%     F_FillArea(mean(Inhibited, 1), sem_, [76, 148, 199]./255, ...
%         1:length(Inhibited))
%     plot(mean(Inhibited, 1), "Color", [76, 148, 199]./255)
% 
%     ylim([0, max(yl)])
%     xlim([1, length(Active)])
%     set(gcf, 'Position', [400, 100, 1200, 500])
%     box off
%     xlabel("Time (Frames)")
%     ylabel("Mean Fluorescence")
%     legend([join(TunningEpochs, ' & '), "", "Excited", "", "Inhibited"])
%     hold off
%     savename = strcat(save_path, "\Tunning - Group mean.pdf");
%     exportgraphics(gcf, savename, "Resolution", 1200)
% 
% % Third figure - Individual-neuron level visualisation
%     close all
%     % Selecting the top five neurons;
%     Sorted = sort(Observed);
%     TopInactive = AllTraces(Observed < Sorted(6), :);
% 
%     % Getting the axes
%     area(double(BinarisedTask).*13, 'EdgeColor','none', ...
%         "FaceAlpha", .3, 'FaceColor', 'k')
%     hold on
%     for i = 1:5
%         plot(TopInactive(i, :) + i, "Color", [76, 148, 199]./255)
%         hold on
%     end
%     Sorted = flip(Sorted);
%     TopActive = AllTraces(Observed > Sorted(6), :);
%     for i = 7:11
% 
%         plot(TopActive(i-6, :)+i, "Color", [212, 100, 66]./255)
%         hold on
%     end
%     hold off
%     ylim([0, 13])
%     yticks([3, 9])
%     yticklabels({"Inhibited", "Excited"}); %#ok<CLARRSTR> 
%     ytickangle(90)
%     xlabel("Time (Frames)")
%     ylabel("Filtered fluorescence")
%     xlim([1, size(AllTraces, 2)])
%     box off
%     Fig_ = gcf;
%     Fig_.Position = [400, 100, 600, 500];
%     savename = strcat(save_path, "\Tunning - Sample Neurons.pdf");
%     exportgraphics(gcf, savename, "ContentType", "vector")
% 
% % Fourth figure - Ratios
%     % Identifying the percentage of neuron tunning
%     close all
%     Properties = [sum(Distances > 1.96, 'All'), ...
%         sum(Distances < -1.96, 'All')];
%     Properties(end+1) = length(Distances) - sum(Properties);
%     pie(Properties)
%     ax = gca();
%     ax.Colormap = [212, 100, 66; 76, 148, 199; 200, 200, 200]./255; 
%     legend(["Excited", "Inhibited", "Irresponsive"], "Location", ...
%         "northeast")
%     F_ = gcf;
%     F_.Position = [400, 100, 390, 320];
%     savename = strcat(save_path, "\Tunning - Ratios.pdf");
%     exportgraphics(gcf, savename, "ContentType", "vector")
% 
% % Fifth figure - Distinction between conditions
%     % Generating the first subplot
%     close all
%     yls = [];
%     subplot(1, 2, 1)
% 
%     boxplot(Output.(ExcitedP_ColName), Output.Sex, ...
%         "Colors", 'k');
%     colours = ...
%         flip(cell2mat(Experiment.Project.Palette(unique(Output.Sex, ...
%         'stable'))), 1);
%     f = findobj(gca, 'Tag', 'Box');
%     for i = 1:length(f)
%         patch(get(f(i),'XData'),get(f(i),'YData'),colours(i,:), ...
%             'FaceAlpha',.6);
%     end
%     h=findobj('LineStyle','--'); set(h, 'LineStyle','-');
%     box off
%     hold on
%     ylabel("Tunned neurons (%)")
%     xlabel("Excited")
%     yls = [yls; ylim()];
% 
%     subplot(1, 2, 2)
%     boxplot(Output.(InhibitedP_ColName), Output.Sex, ...
%         "Colors", 'k');
%     f = findobj(gca, 'Tag', 'Box');
%     for i = 1:length(f)
%         patch(get(f(i),'XData'),get(f(i),'YData'),colours(i,:), ...
%             'FaceAlpha',.6);
%     end
%     h=findobj('LineStyle','--'); set(h, 'LineStyle','-');
%     xlabel("Inhibited")
%     h = gca;
%     h.YAxis.Visible = 'off';
%     box off
%     yls = [yls; ylim()];
%     ylim([min(yls(:, 1)), max(yls(:, 2))])
% 
% 
%     savename = strcat(save_path, "\Tunning - Group Differences.pdf");
%     exportgraphics(gcf, savename, "ContentType", "vector")
%     close all
% 
% %% STEP 7 - CLOSING
% % Saving the report
% savename = strcat(save_path, "\Report.txt");
% writelines(Report,savename);
% % And the data
% savename = strcat(save_path, "\TunningOutput.mat");
% Tunning = Output;
% save(savename, "Tunning");
% 
% % Generating the frequency table
% Responses = table();
% size(SexPerNeuron)
% size(ResponseType)
% size(AnimalPerNeuron)
% Responses.Sex = SexPerNeuron.';
% Responses.Animal = AnimalPerNeuron.';
% Responses.ResponseType = ResponseType.';
% 
% % And saving it
% savename = strcat(save_path, "\SingleNeuronResponses.mat");
% save(savename, "Responses");
% savename_chi = strcat(save_path, "\SingleNeuronResponses.csv");
% writetable(Responses, savename_chi);
% 
% %% STEP 8 - STATS
% % Save Table of Data
% savename = strcat(save_path, "\TunningOutput.csv");
% StatsTable = table();
% StatsTable.Sex = Output.Sex;
% StatsTable.Animal = Output.Animal;
% StatsTable.(ExcitedP_ColName) = Output.(ExcitedP_ColName);
% StatsTable.(InhibitedP_ColName) = Output.(InhibitedP_ColName);
% writetable(StatsTable, savename);
% 
% % Make Table of Variables 
% GroupBy = "Sex";
% Variable1 = ExcitedP_ColName;
% Variable2 = InhibitedP_ColName;
% Id = "Animal";
% DatasetPath = string(savename);
% NewFolderName = "StatResults";
% NewFolderPath = strcat(save_path, "\", NewFolderName);
% DatasetPath_ChiTest = string(savename_chi)
% 
% F_MakeCSVForStat(GroupBy, Variable1, Variable2, Id, DatasetPath, NewFolderPath, DatasetPath_ChiTest);
% 
% % Run R Script of Analysis
% ScriptPath = strcat('', pwd, '\BoxTestAndANOVA.R')
% F_RunRScript(ScriptPath)
% 
% 
% end
% 
% 
