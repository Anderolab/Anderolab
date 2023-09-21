function hb = F_BarPlusError(Means, Error, GroupNames, EpochNames, ...
    ColorDict)

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