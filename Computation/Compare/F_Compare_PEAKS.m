function [Out] = F_Compare_PEAKS(Data1_peak, Data2_peak, RF, GroupBy, ...
    AnimalReg, EpochNames, Visualise, Y_Label, ColorDict,save_path,comparisons,comparisons_colum_name,AnimalsOut_colum)

    
%%
Groups = [];
Animals = [];
Peaks=[];
Animals = unique([AnimalReg{:}]);
Groups = string(GroupBy(Animals).');
Epochs_peak = zeros(length(Animals), 2);
Data_peak={Data1_peak,Data2_peak};

for d_ix = 1:2
    
    % 3 - And through each animal
    for animal = Animals
        
        select_data_peak = Data_peak{d_ix}(AnimalReg{d_ix} == animal, :);
        if size(select_data_peak, 1) ~= 0
            select_data_peak = ...
                select_data_peak(:, ~isnan(select_data_peak(1, :)));
        end

        Epochs_peak(Animals == animal, d_ix) = ...
            sum(sum(select_data_peak==1)) / ...
            ((size(select_data_peak,2)/RF)*size(select_data_peak,1));
    end
end
%end
Animals = Animals.';

% Generating the output table
Out = table(Groups, Animals, Epochs_peak);
% And saving it
savename = strcat(save_path, "\Compare_PEAKS_Table.mat");
save(savename, "Out");
savename_chi = strcat(save_path, "\Compare_PEAKS_Table.csv");
writetable(Out, savename_chi);

%% STATS
% Variables For Grubbs Test
GroupBy = ['Groups'];
Variable1 = ['Epochs_peak_1'];
Variable2 = ['Epochs_peak_2'];
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
NewFolderName = 'StatResults_PEAKS'
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
    
    
    %Read the csv file and identify the outliars
    anova_results = readtable(stats);
    save('anova_results')
    anova_results.Properties.VariableNames{1} = comparisons_colum_name;
    anova_results.Properties.VariableNames{1}
    %if any anova_results.Properties.VariableNames{1}=='SummaryOutliars""'

    [outliar_row, ~] = find(strcmp(anova_results{:,comparisons_colum_name}, 'SummaryOutliars""'))
    % Salta una fila para obtener el primer valor de 'Animal'
    animal_row_start = outliar_row + 1;
    
    % Inicializa un array vacío para almacenar los animales outliers
    outliar_animals = [];
    
    % Sigue buscando hasta que encuentre una fila no válida 
    animal_row = animal_row_start;
    anova_results(animal_row, AnimalsOut_colum)
    
    while animal_row <= size(anova_results, 1) & any(~isnan(str2double(anova_results{animal_row, AnimalsOut_colum}{:})))
        % Extrae el valor del 'Animal' y lo añade al array
        outliar_animals = [outliar_animals; anova_results{animal_row, AnimalsOut_colum}];
    
        % Incrementa el contador de filas
        animal_row = animal_row + 1;
    end
    animal_map = containers.Map('KeyType', 'char', 'ValueType', 'int32');
    for i = 1:length(Animals)
        animal_map(num2str(Animals(i))) = i;
    end
    outliar_animals
    outliar_indices = cellfun(@(x) animal_map(x), num2cell(outliar_animals));
    Grs=unique(Groups);
    Means = [];
    STDs = [];
    N = [];
    for i=1:length(Grs)
        group_index = (Groups == Grs(i));

        % Exclude outlier animals from the group
        valid_animals = setdiff(find(group_index), outliar_indices);
    
        % Compute the mean and std using only the valid animals
        Means = [Means; nanmean(Epochs_peak(valid_animals, :))];
        STDs = [STDs; nanstd(Epochs_peak(valid_animals, :))];
        N = [N; sum(Epochs_peak(valid_animals, :) - ...
            Epochs_peak(valid_animals, :) == 0)];
    end

    SEM = STDs./sqrt(N);

    

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
    
    F_PointsAndMeansLuc(Animals,Epochs_peak, Groups, EpochNames, ColorDict, true,true, 5,4, stats,'PosthocEpochsGroups',comparisons)
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
