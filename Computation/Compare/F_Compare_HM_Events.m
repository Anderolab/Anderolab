function [Out] = F_Compare_HM_Events(Data1, Data2, LowBound, HighBound, ...
    RF, GroupBy, AnimalReg, EpochNames, Visualise, Y_Label, ColorDict,save_path)

%% Actual skeleton for the function
Groups = [];
Animals = [];
Epochs = [];

% 1 - Defining the animals
Animals = unique([AnimalReg{:}]);
Groups = string(GroupBy(Animals).');

    % Changing output size given the number of animals
    Epochs = zeros(length(Animals), 2);

% 2 - Defining single array that stores the data
Data = {Data1, Data2};

% NOTE - IT WOULD BE GOOD TO ENABLE MORE THAN 2 INPUTS FOR THE FUNCTION
% 2 - Looping through each cathegory
for d_ix = 1:2
    
    % 3 - And through each animal
    for animal = Animals

        % Selecting the dataset
        select_data = Data{d_ix}(AnimalReg{d_ix} == animal, :);
        if size(select_data, 1) ~= 0
            select_data = select_data(:, ~isnan(select_data(1, :)));
        end
        subplot(1, 2, 1)
        % Computing the number of events
        num_events = length(F_FindEvents(select_data, LowBound, ...
            HighBound, true));
        xlim([14600, 15600])
        subplot(1, 2, 2)
        num_events = length(F_FindEvents(select_data, LowBound, ...
            HighBound, true));
        xlim([18600, 19600])
        linkaxes(get(gcf,'children'), 'y')
        Epochs(Animals == animal, d_ix) = ...
            num_events / (sum(~isnan(select_data)/RF));
        
    end
end

% Transposing for the sake of generating the output
Animals = Animals.';

% Generating the output table
Out = table(Groups, Animals, Epochs);
% And saving it
savename = strcat(save_path, "\Compare_HM_Events.mat");
save(savename, "Out");
savename_chi = strcat(save_path, "\Compare_HM_Events.csv");
writetable(Out, savename_chi);

%% STATS
% Variables For Grubbs Test
GroupBy = ['Groups'];
Variable1 = ['Epochs_1'];
Variable2 = ['Epochs_2'];
Id = ['Animals'];
DatasetPath = char(savename_chi);
ScriptPath = strcat('', pwd, '\GrubbsTestOutliersEpoch.R'); % R script path

% Run R Script for Grubbs Test of Outliers
threshold = 'FALSE'; % Needed to run correctly the R code
command = ['"' 'C:/PROGRA~1/R/R-43~1.0/bin/Rscript" "' ScriptPath ...
    '" "' DatasetPath '" "' GroupBy '" "' Variable1 '" "' Variable2 ...
    '" "' threshold '"'];
system(command);

% Outliers Table 
Output = readtable("OutliarsTable.csv");

DataTable2 = readtable(DatasetPath);
DataTable2.Animals = string(DataTable2.Animals);
DataTable2.Groups = string(DataTable2.Groups);

% To Remove Outliers 
if width(Output) ~= 1;
    for x = height(Output):-1:1;
        outlier = string(table2array(Output(x,2)));
        [index,~] = find(DataTable2{:,:}==outlier);
        DataTable2(index,:) = []; % deletion
    end 
end

writetable(DataTable2, "DataGlobalOutliersRemoved.csv");

% Variables For Statistical Analysis
NewFolderName = 'StatResults_HM_Events'
NewFolderPath = strcat(save_path, '\', NewFolderName);
DatasetPath = strcat('', pwd, '\DataGlobalOutliersRemoved.csv');
ScriptPath = strcat('', pwd, '\BoxTestAndANOVA.R'); % R script path
ChiTest = 'FALSE';

% Run R Script of Analysis
threshold = 'FALSE'; % Needed to run correctly the R code
command = ['"' 'C:/PROGRA~1/R/R-43~1.0/bin/Rscript" "' ScriptPath ...
    '" "' DatasetPath '" "' GroupBy '" "' Variable1 '" "' Variable2 ...
    '" "' Id '" "' NewFolderPath '" "' ChiTest...
    '" "' threshold '"'];
system(command);

movefile(DatasetPath, NewFolderPath);

%% CALCULATING THE MEANS AND THE STD FOR EACH GROUP

if Visualise

    %% First visualisation - Bars with errors

    Grs = unique(Groups);
    Means = [];
    STDs = [];
    N = [];
    for gr_ix = 1:length(Grs)
        Means = [Means; nanmean(Epochs(Groups == Grs(gr_ix), :))];
        STDs = [STDs; nanstd(Epochs(Groups == Grs(gr_ix), :))];
        N = [N; sum(Epochs(Groups == Grs(gr_ix), :) - ...
            Epochs(Groups == Grs(gr_ix), :) == 0)];
    end
    
    SEM = STDs./sqrt(N);

    % Significance/Non-Significance
    stats = strcat(save_path, "\", NewFolderName, "\SummaryStatistics.csv");
    opts = detectImportOptions(stats);
    columnNames = opts.VariableNames;
    for i = 1:numel(columnNames)
        opts = setvartype(opts, columnNames{i}, 'char');
    end
    stats_table = readtable(stats, opts);
    col_name = stats_table.Properties.VariableNames{2};
    if isempty(stats_table(strcmp(stats_table.(col_name), 'No significant effects'),:)) == false
        csv_check = false;
    elseif isempty(stats_table(strcmp(stats_table.(col_name), '"No significant effects"'),:)) == false
        csv_check = false;
    else 
        csv_check = true; 
    end
    comparisons = { 'Epochs_peak_2:Male-Epochs_peak_2:Female', 'Epochs_peak_2:Male-Epochs_peak_1:Male', ...
    'Epochs_peak_2:Female-Epochs_peak_1:Female','Epochs_peak_1:Male-Epochs_peak_1:Female'};
    
    % Setting y limit register
    max_y = [];
    subplot(1, 2, 1)
    hb =  F_BarPlusErrorLuc(Means.', SEM.', Grs, EpochNames, ColorDict,csv_check,5,stats,'PosthocEpochsGroups',comparisons);
    ylabel(Y_Label)
    max_y = [max_y, max(ylim(gca))];
    title("With incomplete values")
    hold off
    
    %% Second visualisation - Points
    subplot(1, 2, 2)
    
    F_PointsAndMeans(Epochs, Groups, EpochNames, ColorDict, true)
    max_y = [max_y, max(ylim(gca))];
    title("Without incomplete values")
    subplot(1, 2, 1)
    ylim([0, max(max_y)])
    subplot(1, 2, 2)
    ylim([0, max(max_y)])
end

end

