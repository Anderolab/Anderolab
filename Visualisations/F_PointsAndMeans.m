function F_PointsAndMeans(Data, Groups, EpochNames,...
    ColorDict, VisualiseIncomplete)

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

% Computing the means and std for the mean line plot
for g_ix = 1:size(CathegoryName)

    % Setting the colour
    set_col = ColorDict(CathegoryName(g_ix));

    % Computing means and errors
    Mean = mean(Data(Groups == CathegoryName(g_ix), :));
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

% Setting the x labels
xticks([1, length(EpochNames)])
xticklabels(EpochNames)

% Incorporating the legend
legend(Leg)

hold off

end

