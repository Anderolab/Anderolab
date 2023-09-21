function hb = F_BarPlusErrorLuc(Means, Error, GroupNames, EpochNames, ColorDict,csv_logic,pval_colum,csv_name,comparisons_colum_name,comparisons)

% INPUTS
    % Means - rows are different epochs, columns are different groups. The
        % data are the means for individual groups in individual epochs. 
    % Errors - Error bar values in the same configuration as above.
    % Epoch name - Array bearing the names of the epochs represented by
        % each column of Data
    % ColorDict - Dictionary associating each group (key) to a specific
        % colour.
    % VisualiseIncomplete - Boolean. If true, data for animals where
        % measures are incomplete will be represented. These are not
        % represented in the plotted mean.
    %GroupNames - Array bearing the names of the sexes compared.
    %csv_logic=Boolean. If true, there is a csv to read the p-values of the
        %tests performed.
    %pval_colum -Colum in which the p-values are in the input csv.
    %csv_name - Name of the input csv.
    %comparisons_colum_name - Name of the cell just above the comparisons
        %cells. The name has to be without quotation marks.
    %Name of tghe csv cells of the comparisons tests that we want to add 
    % to our plots.


if csv_logic==true
    % Cargar resultados de ANOVA desde archivo CSV
    anova_results = readtable(csv_name);
    
    % Eliminar las comillas dobles al final del nombre de la columna
    anova_results.Properties.VariableNames{1} = comparisons_colum_name;
    
    
    selected_results = table();
    for i = 1:numel(comparisons)
        comparison = comparisons{i}
        [~, index] = ismember(comparison, ...
            table2array(anova_results(:, comparisons_colum_name)))
        
        pval = anova_results{index, pval_colum};  
            % Obtener el valor de p-value de la columna 5
        selected_result = anova_results(index, {comparisons_colum_name});
        selected_result.pval = pval;
        selected_results = [selected_results; selected_result]
    end
    selected_results
    %signif = selected_results.pval
    % Conversión de los p-valores a notación de asteriscos
    %selected_results.pval = cellfun(@str2double, selected_results.pval);
    selected_results.pval_stars = cell(height(selected_results), 1);
    selected_results.pval_stars(selected_results.pval < 0.001) = {'***'};
    selected_results.pval_stars(selected_results.pval >= 0.001 & ...
        selected_results.pval < 0.01) = {'**'};
    selected_results.pval_stars(selected_results.pval >= 0.01 & ...
        selected_results.pval < 0.05) = {'*'};
    selected_results.pval_stars(selected_results.pval >= 0.05) = {'ns'};
    selected_results
    % Generación de las barras
    hb = bar(Means)
    
    
    hold on
    
    % Añadiendo las barras de error
    for k = 1:size(Means,2)
        
        x = hb(k).XData + hb(k).XOffset;
        errorbar(x, Means(:,k), Error(:,k), 'LineStyle', 'none', ...
            'Color', 'k', 'LineWidth', 1);
    end
    height(selected_results)
    % Añadiendo p-valores en notación de asteriscos
    for i = 1:height(selected_results)
        
        
        
        comparison = selected_results.(comparisons_colum_name){i};
        stars = selected_results.pval_stars{i};
        [epoch_a, epoch_b, sex_a, sex_b,group_str] = ...
            extract_groups_from_comparison(comparison);
        disp(group_str);
        disp(['Epochs: ', epoch_a, ' ', epoch_b]);
        
        EpochNames{1} = '1';
        EpochNames{2} = '2';
        EpochNames = cellfun(@char, EpochNames, 'UniformOutput', false)
        sex_indices=[];
        if strcmp(epoch_a, epoch_b)
            % If the epochs are the same, we can just find one index
            epoch_indices = find(ismember(EpochNames, epoch_a));
            % Repeat it to get two indices
            epoch_indices = repmat(epoch_indices, 1, 2)
        else
            % If they're different, find the two unique indices
            epoch_indices = find(ismember(EpochNames, {epoch_a, epoch_b}))
        end
        if strcmp(sex_a, sex_b)
            % If the sexes are the same, we can just find one index
            sex_indices = find(strcmp(GroupNames, sex_a));
            % Repeat it to get two indices
            sex_indices = repmat(sex_indices, 1, 2)
        else
            % If they're different, find the two unique indices
            for i = 1:length({sex_a, sex_b})
                sex_list = {sex_a, sex_b};
                sex_indices = [sex_indices, find(strcmp(GroupNames, ...
                    sex_list{i}))]
            end
        end
        
        % Check if exactly two indices were found for each
        if ~strcmp(epoch_a, epoch_b)
            assert(length(unique(epoch_indices)) == 2);
        end
        if ~strcmp(sex_a, sex_b)
            assert(length(unique(sex_indices)) == 2);
        end
        
    
        % Comprobar si se encontraron exactamente dos índices para cada uno
        assert(length(epoch_indices) == 2);
        assert(length(sex_indices) == 2);
    
        
    
        % Obtiene las coordenadas y para las barras correspondientes
        y_coord = max(Means)
        disp(size(hb(sex_indices(1)).XData));
        disp(size(hb(sex_indices(2)).XData));
        disp(epoch_indices);
    
        % Obtiene las coordenadas x para las barras correspondientes
        % Debe ajustarse de acuerdo con la posición real de las barras en 
        % el gráfico
        size(hb(sex_indices(1)).XData)
        
        size(hb(sex_indices(1)).XOffset)
        
        x_coord = [hb(sex_indices(1)).XData(epoch_indices(1)) + ...
            hb(sex_indices(1)).XOffset, ...
            hb(sex_indices(2)).XData(epoch_indices(2)) + ...
            hb(sex_indices(2)).XOffset];
        
        min(x_coord)
        abs(x_coord(1)-x_coord(2))
        
        y_limits = ylim;
        ylim([y_limits(1),y_limits(2)+(y_limits(2)-y_limits(1))/8]);
        y_limits = ylim; % get the y limits
        
       
        if sex_indices(1)==sex_indices(2) & sex_indices(1)==1
            line(x_coord, [max(y_coord) + 4*(y_limits(2)-y_limits(1))/20 , max(y_coord) + 4*(y_limits(2)-y_limits(1))/20], 'Color', 'black');
            text(mean(x_coord), max(y_coord)+5*(y_limits(2)-y_limits(1))/20, stars, 'HorizontalAlignment', 'center');
        elseif sex_indices(1)==sex_indices(2) & sex_indices(1)==2
            line(x_coord, [max(y_coord) + 8*(y_limits(2)-y_limits(1))/20, max(y_coord) + 8*(y_limits(2)-y_limits(1))/20], 'Color', 'black');
            text(mean(x_coord), max(y_coord)+9*(y_limits(2)-y_limits(1))/20, stars,'HorizontalAlignment', 'center');
        else
            line(x_coord, [max(y_coord) + (y_limits(2)-y_limits(1))/20, max(y_coord) + (y_limits(2)-y_limits(1))/20], 'Color', 'black');
            text(mean(x_coord), max(y_coord)+2*(y_limits(2)-y_limits(1))/20, stars,'HorizontalAlignment', 'center');

        end
        
    end
    
    % Finalizando la figura
    legend(GroupNames)
    set(gca,'xticklabel', EpochNames);
    
    % Pintando
    for b_ix = 1:size(Means, 1)
        col = ColorDict(GroupNames(b_ix));
        hb(b_ix).FaceColor = [col{:}];
    end
    
    hold off

else
        % Generating the bars
    hb = bar(Means);
    
    hold on
    
    % Adding the error bars
    for k = 1:size(Means,2)
    
        x = hb(k).XData + hb(k).XOffset;
        
        % Generating the error bars
        errorbar(x, Means(:,k), Error(:,k), 'LineStyle', 'none', ... 
            'Color', 'k', 'LineWidth', 1);
    end
    
    % Finalising the figure
    legend(GroupNames)
    set(gca,'xticklabel', EpochNames);
    
    % Painting
    for b_ix = 1:size(Means, 1)
        col = ColorDict(GroupNames(b_ix));
        hb(b_ix).FaceColor = [col{:}];
    end
    
    hold off
end

function [epoch_a, epoch_b, sex_a, sex_b,group_str] = ...
        extract_groups_from_comparison(comparison)
    % Split the comparison string into parts
    group_str = split(comparison, '-');
    group_a = split(group_str{1}, ':');
    group_b = split(group_str{2}, ':');
    
    % Extract the epochs and sexes
    epoch_a = strtrim(group_a{1});
    epoch_b = strtrim(group_b{1});
    sex_a = strtrim(group_a{2});
    sex_b = strtrim(group_b{2});
    
   % Check if there's '1' or '2' in epoch_a and epoch_b
    if contains(epoch_a, '_1')
        epoch_a = '1';
    elseif contains(epoch_a, '_2')
        epoch_a = '2';
    end

    if contains(epoch_b, '_1')
        epoch_b = '1';
    elseif contains(epoch_b, '_2')
        epoch_b = '2';
    end
end
end