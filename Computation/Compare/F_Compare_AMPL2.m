function [Out] = F_Compare_AMPL2(Data1, Data2, GroupBy, ...
    AnimalReg, EpochNames, Visualise, Y_Label, ColorDict, PeakType, save_path)
%%
Groups = [];
Animals = [];
Peaks=[];
Animals = unique([AnimalReg{:}]);
Groups = string(GroupBy(Animals).');
Epochs_peak = zeros(length(Animals), 2);
Data={Data1,Data2};

for d_ix = 1:2
    
    % 3 - And through each animal
    for animal = Animals
        
        select_data = Data{d_ix}(AnimalReg{d_ix} == animal, :);
        if PeakType == "Raw"

            peaks = islocalmax(select_data, 2);

        else
            peaks = F_ExtractFrequency(select_data, 30, 30);
            peaks = peaks.peaks;
           
            peaks(isnan(select_data)) = NaN;
        end
        
        amps = peaks.*select_data;

        % Calcule the number of ones
        ones_per_row = nansum(peaks, 2);
       
        
        % Calcular la media de cada fila
        row_sum = nansum(amps, 2);
        
        % Calcular la suma de las medias de las filas, dividir entre por el número de 1s y luego por el número de filas
        Epochs(Animals == animal, d_ix) = ...
            nansum(row_sum ./ ones_per_row) / size(select_data, 1);
    end
end
%end
Animals = Animals.';

% Generating the output table
Out = table(Groups, Animals, Epochs)
% And saving it
savename = strcat(save_path, "\Compare_AMPL_Table.mat");
save(savename, "Out");
savename_chi = strcat(save_path, "\Compare_AMPL_Table.csv");
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
NewFolderName = 'StatResults_AMPL'
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
    Grs=unique(Groups);
    Means = [];
    STDs = [];
    N = [];
    for i=1:length(Grs)
        Means=[Means;nanmean(Epochs(Groups == Grs(i), :))];
        STDs = [STDs; nanstd(Epochs(Groups == Grs(i), :))];
        N = [N; sum(Epochs(Groups == Grs(i), :) - ...
            Epochs(Groups == Grs(i), :) == 0)];
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
    comparisons = {'Epochs_2:Male-Epochs_2:Female', 'Epochs_2:Male-Epochs_1:Male', ...
    'Epochs_2:Female-Epochs_1:Female','Epochs_1:Male-Epochs_1:Female'};

    max_y = [];
    subplot(1, 2, 1)
    hb =  F_BarPlusErrorLuc(Means.', SEM.', Grs, EpochNames, ColorDict,csv_check,5,stats,'PosthocEpochsGroups',comparisons);
    ylabel('Frequency peaks')
    max_y = [max_y, max(ylim(gca))];
    title("With incomplete values")
    hold off

    %% Second visualisation - Points
    subplot(1, 2, 2)
    
    F_PointsAndMeans(Epochs, Groups, EpochNames, ColorDict, true)
    max_y = [max_y, max(ylim(gca))];
    title("Without incomplete values")
    
    subplot(1, 2, 1)
    ylabel(Y_Label)
    ylim([0, max(max_y)])
    subplot(1, 2, 2)
    ylim([0, max(max_y)])
    f = gcf;
    f.Position = [100 100 540 400];
end

end