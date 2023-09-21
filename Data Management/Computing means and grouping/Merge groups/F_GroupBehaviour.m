function [FreezData, Means, SDs, N] = F_GroupBehaviour(Experiment, ...
    FreezTable, EpochsOfInterest, GroupingVariable, Session);

% Slicing
FreezData = FreezTable(:, ...
    contains(string(FreezTable.Properties.VariableNames), ...
    EpochsOfInterest));
FreezData.Animal = FreezTable.Animal;

%% Grouping

% Determining group of each animal
    Groups = string(Experiment.Project.(GroupingVariable));

    a_ix = [];
    for a = FreezData.Animal.'
        a_ix = [a_ix; sscanf(a,'M%d')];
    end
    FreezData.(GroupingVariable) = Groups(a_ix).';

% Computing group mean
    
    
    Groups = unique(Groups);

    % Empty storage
    Means = zeros(length(Groups), length(EpochsOfInterest));
    SDs = Means;
    N = zeros(1, length(Groups));
    
    
    c = 1; % Counter
    for gr = Groups

        % Extracting freex per group
        gr_freez = table2array(FreezData(FreezData.(GroupingVariable) == gr, ...
            1:length(EpochsOfInterest)));

        % Computing
        Means(c, :) = mean(gr_freez, 1);
        SDs(c, :) = std(gr_freez, [], 1);
        N(c) = size(gr_freez, 1);
        c = c+1;
    end

%% Stats - (FreezData)

FreezData = [FreezData(:,end) FreezData(:,end-1) FreezData(:,2:end-2)];
writetable(FreezData, "FreezData.csv");
F_ToneOutlierRemoval('FreezData.csv', 'Sexes', 'Tone', 'Animal');
F_GlobalOutlierRemovalPlusSTATS('FreezData.csv', 'Sexes', 'Tone', 'Animal');

new_folder = strcat('StatsResults', Session);
mkdir(new_folder);
movefile("StatResultsGlobalOutliers", new_folder);
movefile("StatResultsOutliersByTone", new_folder);
movefile("OutliarsTable.csv", new_folder);

%% Visualising (optional)
for i = 1:(c-1)
    col = Experiment.Project.Palette(Groups(i));
    F_FillArea(Means(i, :), SDs(i, :)./sqrt(N(i)), col{:}, 1:length(Means(i, :)))
    hold on
    plot(Means(i, :), "Color", col{:}, "LineWidth", 2)
    
    hold on
end
xticks(1:length(EpochsOfInterest))
xticklabels(regexprep(EpochsOfInterest, '_', ' '))
end

