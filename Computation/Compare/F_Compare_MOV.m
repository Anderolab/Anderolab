function [Out] = F_Compare_MOV(Data1, Data2, Sexes, AnimalReg, ...
    EpochNames, Visualise, Y_Label, ColorDict,LowBound,HighBound,save_path)
%% Counting
Counts = [];
Animals = [];
Epochs = [];

% Defining animals
Animals = unique([AnimalReg{:}]);

Groups = string(Sexes(Animals)).';

% Defining the output
Epochs = zeros(length(Animals), 2);

% Generating single array for storage
Data = {Data1, Data2};

for d_ix = 1:2

    for animal = Animals
        select_data = Data{d_ix}(AnimalReg{d_ix} == animal, :);
        % Computing
        high_mov = sum(select_data > LowBound & select_data < HighBound);
        non_nans = sum(~(isnan(select_data)));
        ratio = high_mov/non_nans;

        % Saving
        Epochs(Animals == animal, d_ix) = ratio;
    end
end

% Transposing
Animals = Animals.';

% Making the Table 
Out = table(Groups, Animals, Epochs);
% And saving it
savename = strcat(save_path, "\Compare_MOV_Table.mat");
save(savename, "Out");
savename_chi = strcat(save_path, "\Compare_MOV_Table.csv");
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
NewFolderName = 'StatResults_MOV'
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
    comparisons = { 'Epochs_2:Male-Epochs_2:Female', 'Epochs_2:Male-Epochs_1:Male', ...
    'Epochs_2:Female-Epochs_1:Female','Epochs_1:Male-Epochs_1:Female'};

    % Setting y limit register
    max_y = [];
    min_y = [];
    subplot(1, 2, 1)
    hb = F_BarPlusErrorLuc(Means.', SEM.', Grs, EpochNames, ColorDict,csv_check,5,stats,'PosthocEpochsGroups',comparisons);
    ylabel(Y_Label)
    max_y = [max_y, max(ylim(gca))];
    min_y = [min_y, min(ylim(gca))];
    title("With incomplete values")
    hold off
    
    
    %% Second visualisation - Points
    subplot(1, 2, 2)
    
    F_PointsAndMeans(Epochs, Groups, EpochNames, ColorDict, true)
    max_y = [max_y, max(ylim(gca))];
    min_y = [min_y, min(ylim(gca))];
    title("Without incomplete values")
    subplot(1, 2, 1)
    ylim([max(min_y), max(max_y)])
    subplot(1, 2, 2)
    ylim([max(min_y), max(max_y)])
end

end

