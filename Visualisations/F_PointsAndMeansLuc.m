function F_PointsAndMeansLuc(Animals,Data, Groups, EpochNames, ColorDict, VisualiseIncomplete, csv_logic, pval_colum,AnimalsOut_colum, csv_name, comparisons_colum_name, comparisons)

% INPUTS
    % Data - rows are individual animals, columns are the different
        % epochs.
    % Groups - List of string containing the cathegory name to which each
        % animal belongs to.
    % Epoch name - Array bearing the names of the epochs represented by
        % each column of Data
    % ColorDict - Dictionary associating each group (key) to a specific
        % colour.
    % VisualiseIncomplete - Boolean. If true, data for animals where
        % measures are incomplete will be represented. These are not
        % represented in the plotted mean.

if csv_logic~=true
    % Creating an empty legend to populate
    Leg = [];
    
    CathegoryName = unique(Groups);
    
    % Identifyng the nans
    NonNaN_IX = sum((Data-Data) == 0, 2) == size(Data, 2);
    
    % Visualising the non-complete values
    if VisualiseIncomplete == true
        
        % Finding indexes bearing NaNs and selicing accordingly
        NaN_IX = sum((Data-Data) == 0, 2) < size(Data, 2);
        NaN_Dat = Data(NaN_IX, :);
        NaN_Cath = Groups(NaN_IX, :);
        
        % Plotting empty values
        for n_ix = 1:length(NaN_Cath)
    
            % Cathegory-specific colour selection
            col = ColorDict(NaN_Cath(n_ix));
    
            % Identifying epoch bearing the non-nan values
            x = find((NaN_Dat(n_ix, :)-NaN_Dat(n_ix, :) == 0) == 1);
    
            % Visualising
            scatter(x, NaN_Dat(n_ix, x), 70, col{:}, "Marker", '+');
    
            % Updating the legend
            Leg = [Leg, ""];
            hold on
        end
    end
    
    % Slicing the complete data
    Groups = Groups(NonNaN_IX, :);
    Data = Data(NonNaN_IX, :);
    
    % Plotting
    for p_ix = 1:sum(NonNaN_IX)
    
        % Cath-specific colour selection
        set_col = ColorDict(Groups(p_ix));
    
        % Plotting points and lines
        plot(Data(p_ix, :), "Color", set_col{:}, "LineStyle", ":")
        hold on
        scatter(1:length(Data(p_ix, :)), Data(p_ix, :), [], set_col{:})
        Leg = [Leg, "", ""];
    end
    
    % Setting the limits
    xlim([0.75, length(CathegoryName)+.25])
    hold on
    % Initialize the matrix
    Means = nan(2, size(Data, 2)); 
    % Computing the means and std for the mean line plot
    for g_ix = 1:size(CathegoryName)
    
        % Setting the colour
        set_col = ColorDict(CathegoryName(g_ix));
    
        % Computing means and errors
        Mean = mean(Data(Groups == CathegoryName(g_ix), :))
        Means(g_ix,:) = Mean;
        STDs = std(Data(Groups == CathegoryName(g_ix), :));
        SEM = STDs./sqrt(size(Data(Groups == CathegoryName(g_ix), :), 1));
    
        % Plotting
        errorbar(1:length(Data(p_ix, :)), Mean, SEM, 'LineStyle', 'none', ... 
            'Color', 'w', 'LineWidth', 5);
        errorbar(1:length(Data(p_ix, :)), Mean, SEM, 'LineStyle', 'none', ... 
            'Color', set_col{:}, 'LineWidth', 1.5);
        hold on
        plot(Mean, "Color", 'w', "LineWidth", 4)
        plot(Mean, "Color", set_col{:}, "LineWidth", 3)
        scatter(1:length(Data(p_ix, :)), Mean, 90, set_col{:}, "filled")
        Leg = [Leg, "", "", "", CathegoryName(g_ix), ""];
    
    
    end
end
hold on 
CathegoryName = unique(Groups);
if csv_logic==true
    % Cargar resultados de ANOVA desde archivo CSV
    anova_results = readtable(csv_name);
    % Eliminar las comillas dobles al final del nombre de la columna
    anova_results.Properties.VariableNames{1} = comparisons_colum_name;
     % Busca la fila que contiene 'SummaryOutliars'
    [outliar_row, ~] = find(strcmp(anova_results{:,comparisons_colum_name}, 'SummaryOutliars""'))
    
    % Salta una fila para obtener el primer valor de 'Animal'
    animal_row_start = outliar_row + 1;
    
    % Inicializa un array vacío para almacenar los animales outliers
    outliar_animals = [];
    
    % Sigue buscando hasta que encuentres una fila no válida (por ejemplo, una fila en blanco)
    animal_row = animal_row_start;
    anova_results(animal_row, AnimalsOut_colum)
    while animal_row <= size(anova_results, 1) & any(~isnan(str2double(anova_results{animal_row, AnimalsOut_colum}{:})))
        % Extrae el valor del 'Animal' y lo añade al array
        outliar_animals = [outliar_animals; anova_results{animal_row, AnimalsOut_colum}];
    
        % Incrementa el contador de filas
        animal_row = animal_row + 1;
    end
    outliar_animals
    % Creating an empty legend to populate
    Leg = [];
    
    
    
    % Identifyng the nans
    NonNaN_IX = sum((Data-Data) == 0, 2) == size(Data, 2);
    
    % Visualising the non-complete values
    if VisualiseIncomplete == true
        
        % Finding indexes bearing NaNs and selicing accordingly
        NaN_IX = sum((Data-Data) == 0, 2) < size(Data, 2);
        NaN_Dat = Data(NaN_IX, :);
        NaN_Cath = Groups(NaN_IX, :);
        
        % Plotting empty values
        for n_ix = 1:length(NaN_Cath)
    
            % Cathegory-specific colour selection
            col = ColorDict(NaN_Cath(n_ix));
    
            % Identifying epoch bearing the non-nan values
            x = find((NaN_Dat(n_ix, :)-NaN_Dat(n_ix, :) == 0) == 1);
    
            % Visualising
            scatter(x, NaN_Dat(n_ix, x), 70, col{:}, "Marker", '+');
    
            % Updating the legend
            Leg = [Leg, ""];
            hold on
        end
    end
    
    % Slicing the complete data
    Groups = Groups(NonNaN_IX, :);
    Data = Data(NonNaN_IX, :);
    Non=NonNaN_IX
    outliar_animals_num = str2double(outliar_animals);  % Convertir a números

    % Plotting
    for p_ix = 1:sum(NonNaN_IX)

        if ismember(Animals(p_ix), outliar_animals_num)
             continue;  % Skip this iteration if it is an outliar
        end
    
        % Cath-specific colour selection
        set_col = ColorDict(Groups(p_ix));
    
        % Plotting points and lines
        plot(Data(p_ix, :), "Color", set_col{:}, "LineStyle", ":")
        hold on
        scatter(1:length(Data(p_ix, :)), Data(p_ix, :), [], set_col{:})
        Leg = [Leg, "", ""];
    end
    
    % Setting the limits
    xlim([0.75, length(CathegoryName)+.25])
    hold on
    % Initialize the matrix
    Means = nan(2, size(Data, 2)); 
    % Computing the means and std for the mean line plot
    animal_map = containers.Map('KeyType', 'char', 'ValueType', 'int32');
    for i = 1:length(Animals)
        animal_map(num2str(Animals(i))) = i;
    end
    outliar_indices = cellfun(@(x) animal_map(x), num2cell(outliar_animals));
    for g_ix = 1:size(CathegoryName)
    
        % Setting the colour
        set_col = ColorDict(CathegoryName(g_ix));
        % Create an index of the animals for this category
        category_index = (Groups == CathegoryName(g_ix));
    
        % Exclude outlier animals from the category
        valid_animals = setdiff(find(category_index), outliar_indices);
        
        % Compute the mean and std using only the valid animals
        Mean = mean(Data(valid_animals, :));
        Means(g_ix,:) = Mean;
        STDs = std(Data(valid_animals, :));
        SEM = STDs./sqrt(length(valid_animals));
    
        % Plotting
        errorbar(1:length(Data(p_ix, :)), Mean, SEM, 'LineStyle', 'none', ... 
            'Color', 'w', 'LineWidth', 5);
        errorbar(1:length(Data(p_ix, :)), Mean, SEM, 'LineStyle', 'none', ... 
            'Color', set_col{:}, 'LineWidth', 1.5);
        hold on
        plot(Mean, "Color", 'w', "LineWidth", 4)
        plot(Mean, "Color", set_col{:}, "LineWidth", 3)
        scatter(1:length(Data(p_ix, :)), Mean, 90, set_col{:}, "filled")
        Leg = [Leg, "", "", "", CathegoryName(g_ix), ""];
    
    end
    
    
    
    selected_results = table();
    for i = 1:numel(comparisons)
        comparison = comparisons{i};
        [~, index] = ismember(comparison, table2array(anova_results(:, comparisons_colum_name)));
        
        pval = anova_results{index, pval_colum};  % Obtener el valor de p-value de la columna 5
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
    selected_results.pval_stars(selected_results.pval >= 0.001 & selected_results.pval < 0.01) = {'**'};
    selected_results.pval_stars(selected_results.pval >= 0.01 & selected_results.pval < 0.05) = {'*'};
    selected_results.pval_stars(selected_results.pval >= 0.05) = {'ns'};
    selected_results
    % Generación de las barras
    for i = 1:height(selected_results)
        
        
        
        comparison = selected_results.(comparisons_colum_name){i};
        stars = selected_results.pval_stars{i};
        [epoch_a, epoch_b, sex_a, sex_b,group_str] = extract_groups_from_comparison(comparison);
        disp(group_str);
        disp(['Epochs: ', epoch_a, ' ', epoch_b]);
        
        
        EpochNames
        epoch_a
        epoch_b
        
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
            sex_indices = find(strcmp(CathegoryName, sex_a));
            % Repeat it to get two indices
            sex_indices = repmat(sex_indices, 1, 2)
        else
            % If they're different, find the two unique indices
            for i = 1:length({sex_a, sex_b})
                sex_list = {sex_a, sex_b};
                sex_indices = [sex_indices, find(strcmp(CathegoryName, sex_list{i}))]
                sex_indices
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
        y_samesex=max(max(Data))
        

        % y_coord = max(Means)
        % disp(size(hb(sex_indices(1)).XData));
        % disp(size(hb(sex_indices(2)).XData));
        % disp(epoch_indices);
    
        % Obtiene las coordenadas x para las barras correspondientes
        % Debe ajustarse de acuerdo con la posición real de las barras en el gráfico
        % size(hb(sex_indices(1)).XData)
        % 
        % size(hb(sex_indices(1)).XOffset)
        
        % x_coord = [hb(sex_indices(1)).XData(epoch_indices(1)) + hb(sex_indices(1)).XOffset,hb(sex_indices(2)).XData(epoch_indices(2)) + hb(sex_indices(2)).XOffset]
        % 
        % min(x_coord)
        % abs(x_coord(1)-x_coord(2))
        % 
        y_limits = ylim;
        ylim([y_limits(1),y_limits(2)+(y_limits(2)-y_limits(1))/8]);
        y_limits = ylim; % get the y limits
        x_limits = xlim;
        xlim([x_limits(1)-(x_limits(2)-x_limits(1))/25,x_limits(2)+(x_limits(2)-x_limits(1))/35]);
        x_limits = xlim;
        Means
        Means(sex_indices(1),epoch_indices(1))
        Means(sex_indices(2),epoch_indices(2))
        array_text=[Means(sex_indices(1),epoch_indices(1)),Means(sex_indices(2),epoch_indices(2))]
        group1_x_position = epoch_indices(1);
        group2_x_position = epoch_indices(2);
        x_coord_dif_sex=epoch_indices(1);
        x_coord=[group1_x_position,group2_x_position]
        color_male=ColorDict('Male');
        color_female=ColorDict('Female');
        color_male=color_male{1};
        color_female=color_female{1};
        if sex_indices(1)==sex_indices(2) & sex_indices(1)==1
            line(x_coord, [y_samesex + 4*(y_limits(2)-y_limits(1))/20 , y_samesex + 4*(y_limits(2)-y_limits(1))/20], 'Color', 'black');
            text(mean(x_coord), y_samesex+5*(y_limits(2)-y_limits(1))/20, stars,'Color',color_female, 'HorizontalAlignment', 'center');
        elseif sex_indices(1)==sex_indices(2) & sex_indices(1)==2
            line(x_coord, [y_samesex + 8*(y_limits(2)-y_limits(1))/20, y_samesex + 8*(y_limits(2)-y_limits(1))/20], 'Color', 'black');
            text(mean(x_coord), y_samesex+9*(y_limits(2)-y_limits(1))/20, stars,'Color',color_male,'HorizontalAlignment', 'center');
        elseif epoch_indices==1
            line(epoch_indices-(x_limits(2)-x_limits(1))/19 , [Means(sex_indices(1),epoch_indices(1)), Means(sex_indices(2),epoch_indices(2))], 'Color', 'black');
            text(epoch_indices(1)-(x_limits(2)-x_limits(1))/12, mean(array_text), stars,'HorizontalAlignment', 'center');
        elseif epoch_indices==2
            line(epoch_indices+(x_limits(2)-x_limits(1))/19 , [Means(sex_indices(1),epoch_indices(1)), Means(sex_indices(2),epoch_indices(2))], 'Color', 'black');
            text(epoch_indices(1)+(x_limits(2)-x_limits(1))/12, mean(array_text), stars,'HorizontalAlignment', 'center');
        end
        
        end
    end
   
    hold on

% Setting the x labels
xticks([1, length(EpochNames)])
xticklabels(EpochNames)

% Incorporating the legend
legend(Leg)

hold off
function [epoch_a, epoch_b, sex_a, sex_b,group_str] = extract_groups_from_comparison(comparison)
    % Split the comparison string into parts
    group_str = split(comparison, '-')
    group_a = split(group_str{1}, ':')
    group_b = split(group_str{2}, ':')
    
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

