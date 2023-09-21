%% LOAD THE REQUIRED DATASETS
% 1 - Adding the project path to the searchpath
% Experiment properties
output_path = "C:\Users\1657711\OneDrive - UAB\" + ...
    "Documentos\Código Joaquin\TestEnvironment\TestLoad";

% Adding functions to search path
addpath(genpath(output_path))
load("ExperimentData.mat")


trial = "FC";
fig_filename = strcat(trial, ".pdf");
Field = 'Filt';
Epochs = ["Tone FC1", "Tone FC2", "Tone FC3", "Tone FC4", "Tone FC5"];
% Epochs = ["Tone FE1", "Tone FE2", "Tone FE3", "Tone FE4", "Tone FE5", ...
%     "Tone FE6", "Tone FE7", "Tone FE8", "Tone FE9", "Tone FE10", ...
%     "Tone FE11", "Tone FE12", "Tone FE13", "Tone FE14", "Tone FE15"];
load("Contrast.mat")
load('batlow.mat')
Session_data = F_PreClustering(Experiment, trial, ...
    Epochs, true, -30)




%% SECTION 1 - PREPARING AND PRE-PROCESSING THE DATA

%% Selecting neurons with possible distorting factors
Session = [];
Clustering_data = [];
% LOADING FC
% trial = "FC";
Field = 'Filt';
% Epochs = ["Tone FC1", "Tone FC2", "Tone FC3", "Tone FC4", "Tone FC5"] + ...
%     " " + trial;
% Experiment.(trial).Task.Titles = ...
%     string(Experiment.(trial).Task.Titles)+" "+trial;
Session_data = F_PreClustering(Experiment, trial, ...
    Epochs, true, -30);
Session = [Session, repelem(trial, size(Session_data.Flu, 1))];
Clustering_data = [Clustering_data, Session_data];

%% SECTION 2 - IDENTIFYING THE IDEAL CLUSTER NUMBER
%% Runing multiple kmeans tests to identify best number of ks
BestK = F_GetBestK(Session_data.Flu, 200);                               % NOTE FOR JOAQUÍN - ATM performing all options, rewrite as optimisation function for speed purposes
K_Ixes = 1:BestK;
Colors = repmat(Contrast, ceil(BestK/size(Contrast, 1)), 1);
exportgraphics(gcf, fig_filename, "Append", true, "ContentType", 'vector')
%% SECTION 3 - ACTUALLY PERFORMING KMEANS
[Session_data.Clusters] = kmeans(Session_data.Flu, BestK, 'MaxIter', 500);
F_Table = tabulate(Session_data.Clusters);
[Percentage, Ix] = sort(F_Table(:, 3), 'descend')
ColDict = dictionary(K_Ixes(Ix), K_Ixes);
ColorData = Colors(ColDict(Session_data.Clusters), :);

%% Assessing within cluster variability
[Var, ~, Centroids, Klu, N] = F_VarPerK(Session_data.Flu, ...
    Session_data.Clusters);
[~, ix_] = sort(Klu);


% Viewing the variance per cluster
F_ViewDistribution(Var, "Nº clusters", "Within cluster \sigma^2")
exportgraphics(gcf, fig_filename, "Append", true, "ContentType", 'vector')
% Viewing the n per each cluster
F_ViewDistribution(N, "Nº clusters", "N per cluster")
exportgraphics(gcf, fig_filename, "Append", true, "ContentType", 'vector')
%% Visualising the data in the UMAP plane
Clustering_data=Session_data
umap = run_umap(gather(Session_data.Flu), 'n_neighbors', 10, 'min_dist', .5); %, 
close all
%

colormap(batlow)

subplot(1, 26, [1:10])
F_HeatScatter(umap(:, 1), umap(:, 2), 40)
xl = xlim();
yl = ylim();
colorbar
xlabel("UMAP dimension 1")
ylabel("UMAP dimension 2")

subplot(1, 26, [1:9]+12)
scatter(umap(:, 1), umap(:, 2), 7, ColorData, "filled")
xlabel("UMAP dimension 1")
ylabel("UMAP dimension 2")
xlim(xl)
ylim(yl)


% Selecting most common types
subplot(1, 26, [24, 26])
CommonResp = Centroids(Ix(1:12), :);
ys = [1:12]./3
ylabels = []
for kluster = 1:12
    k  = Ix(kluster);
    ylabels = [ylabels, k];
    y = mean(Clustering_data.Flu(Clustering_data.Clusters == k, :), 1);
    y = (y-median(y));
    plot([1:size(Clustering_data.Flu, 2)]./30, y + kluster/3, ...
        'Color', Colors(ColDict(k), :), 'LineWidth', 2);
    hold on
end
hold off
yticks(ys)
xlabel('Time (s)')
ylabel('Cluster \DeltaF')
yticklabels(string(ylabels))
box off
f_ = gcf;
f_.Position = [1033, 457, 1090, 360];
exportgraphics(f_, "UMAP view of clusters_GLOBAL.pdf", 'ContentType', ...
    'vector')
exportgraphics(gcf, fig_filename, "Append", true, "ContentType", 'vector')
%% SECTION 4 - PERFORMING HERARCHICAL CLUSTERING

% Normalising the centroids
NormCents = ...
    (Centroids-min(Centroids, [], 2))./ ...
    (max(Centroids, [], 2) -min(Centroids, [], 2));
Distances = pdist(NormCents); % Computing distances between the centroids
Tree = linkage(Distances, "ward"); % And generating a tree given proximity
cutoff = median([Tree(end-1,3) Tree(end,3)]) % Setting the cutoff

% Generating clustering given the herarchy tree
H_Clusters = cluster(Tree, 'maxclust', 2);
leafOrder = optimalleaforder(Tree, Distances);
cutoff = median([Tree(end-1,3) Tree(end,3)])

%% Generating output for chi square
Sex = string(Experiment.Project.Sexes);
Clustering_data.Sex = Sex(double(Clustering_data.Animal)).';
Clustering_data.Supercluster = H_Clusters(Clustering_data.Clusters);
writetable(Clustering_data, "Clustering_GLOBAL.csv")


%% Extracting relative frequencies

subplot(7, 12, [1, 2, 3, 4, [11, 12, 13, 14]+2, [21, 22, 23, 24]+4, ...
    [31, 32, 33, 34]+6])
    hLines = dendrogram(Tree, 'ColorThreshold', cutoff, ...
        "Reorder", leafOrder);
    hold on
    yline(cutoff, 'LineWidth', 1.2, 'Color', 'r', 'LineStyle', ':')
    ylabel('Ward centroid-centroid distance')
    xlabel('Cluster')
    hold off

% Associating with leaf order
cluster_order = unique(H_Clusters(leafOrder), 'stable');

% To maintain equal axes
yls = [];
subplot(7, 12, [41, 42, 53, 54] + 20)
    % Mean supercluster fluorescence
    supercluster_f = mean(NormCents(H_Clusters == cluster_order(1), :));
    area([1:length(supercluster_f)]./30, supercluster_f)
    yls = [yls, ylim()]
    xlabel("Time")
    ylabel('Supercluster centroid')


subplot(7, 12, [43, 44, 55, 56] + 20)
    supercluster_f = mean(NormCents(H_Clusters == cluster_order(2), :));
    area([1:length(supercluster_f)]./30, supercluster_f)
    yticks([])
    yls = [yls, ylim()]
    xlabel("Time")

subplot(7, 12, [43, 44, 55, 56] + 20)
    ylim([min(yls, [], 'all'), max(yls, [], 'all')])
subplot(7, 12, [41, 42, 53, 54] + 20)
    ylim([min(yls, [], 'all'), max(yls, [], 'all')])
%
subplot(7, 12, [[6:8],  [6:8]+12, [6:8]+24])

    Frequencies = F_FrequencyTable(Clustering_data.Epoch,...
        Clustering_data.Clusters, true);
    
    % Generating the table
    bg_ = bar(table2array(Frequencies(:, 2:end)).', 'stacked')
    Frequencies.Properties.VariableNames(2:end)
    % Labelling it
    xticks([1:length(Frequencies.Properties.VariableNames(2:end))])
    xticklabels(Frequencies.Properties.VariableNames(2:end))
    xlim([.5, size(Frequencies, 2)-1+.5])
    ylim([0, 100])

    % Changing colors
    for k = 1:BestK
        bg_(k).FaceColor = Colors(k, :);
    end
    
    ylabel("% of cases")
    

subplot(7, 12, [[6:8]+12*4, [6:8]+12*5, [6:8]+12*6])

    % Generating the table
    Frequencies = F_FrequencyTable(Clustering_data.Sex,...
        Clustering_data.Clusters, true);
    bg_ = bar(table2array(Frequencies(:, 2:end)).', 'stacked')

    % Labelling it
    xticklabels(Frequencies.Properties.VariableNames(2:end))
    xlim([.5, size(Frequencies, 2)-1+.5])
    ylim([0, 100])
    
    % Changing colors
    for k = 1:BestK
        bg_(k).FaceColor = Colors(k, :);
    end
    ylabel("% of cases")


subplot(7, 12, [[10:12],  [10:12]+12, [10:12]+24])
    Frequencies = F_FrequencyTable(Clustering_data.Epoch,...
        Clustering_data.Supercluster, true)
    
    % Generating the table
    bar(table2array(Frequencies(:, 2:end)).', 'stacked')
    
    % Labelling it
    xticklabels(Frequencies.Properties.VariableNames(2:end))
    xlim([.5, size(Frequencies, 2)-1+.5])
    ylim([0, 100])
    
    ylabel("% of cases")

subplot(7, 12, [[10:12]+12*4, [10:12]+12*5, [10:12]+12*6])
    Frequencies = F_FrequencyTable(Clustering_data.Sex,...
        Clustering_data.Supercluster, true)
    
    % Generating the table
    bar(table2array(Frequencies(:, 2:end)).', 'stacked')
    
    % Labelling it
    xticklabels(Frequencies.Properties.VariableNames(2:end))
    xlim([.5, size(Frequencies, 2)-1+.5])
    ylim([0, 100])
    
    ylabel("% of cases")

f_ = gcf;
f_.Position = [680, 470, 970, 500];

exportgraphics(gcf, fig_filename, "Append", true, "ContentType", 'vector')
%%
close all
% Looking at it on an animal-animal basis
    SexEpoch_F = {};
    Frequency_ = [];
    Sex = [];
    Animal = unique(Session_data.Animal).';
    means = [];
    sds = [];
    ns = [];
    c = 1;
    for gr = unique(Clustering_data.Sex).'
        gr
        gr_data = Clustering_data(Clustering_data.Sex == gr, :);
        gr_f = [];
        for an = unique(gr_data.Animal).'
            an_data = gr_data(gr_data.Animal == an, :);
            Frequency = table2array(F_FrequencyTable(an_data.Epoch,...
                an_data.Supercluster, true));
            gr_f = [gr_f; Frequency(2, 2:end)]
        end
        gr_f
        Frequency_ = [Frequency_; gr_f]
        Sex = [Sex; repelem(gr, size(gr_f, 1)).']
        SexEpoch_F{c} = gr_f
        c = c+1;
        means = [means; mean(gr_f, 1)];
        sds = [sds; std(gr_f, [], 1)];
        ns = [ns, size(gr_f, 1)];
    end
    DataForStatistics = table();
    DataForStatistics.Animal = Animal.';
    DataForStatistics.Sex = Sex
    DataForStatistics.Frequency = Frequency_

%%
    % Hi hunnn was wondering, like, no worries, feel free to do whatevs but 
    % like OMG if you could possibly like do perhaps some stats around here
    % (only if u want, like put them wherever) loke OMG you would like super
    % save me cause lol I'm such an idiot like I can't but ofc no pressure
    % like it's fine if not, i can do it, hun u don't get paid or anything but
    % even if you did, thank u <3

    % ok 

    writetable(DataForStatistics, "DataForStatistics.csv");

    % Variables For Grubbs Test
    GroupBy = ['Sex'];
    Variable1 = ['Frequency'];
    ScriptPath = strcat('', pwd, '\GrubbsTestOutliersMock.R'); % R script path
    DataTable = char("DataForStatistics.csv")
    
    % Run R Script of Analysis
    threshold = 'FALSE'; % Needed to run correctly the R code
    command = ['"' 'C:/PROGRA~1/R/R-43~1.0/bin/Rscript" "' ScriptPath ...
        '" "' DataTable '" "' GroupBy '" "' Variable1 '" "' threshold '"'];
    system(command);
    
    % Outliers Table 
    Output = readtable("OutliarsTable.csv");
    
    DataTable2 = readtable(DataTable);
    DataTable2.Animal = string(DataTable2.Animal);
    DataTable2.Sex = string(DataTable2.Sex);
    
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
    Id = ['Animal'];
    DataGlobalOutliers = strcat('', pwd, '\DataGlobalOutliersRemoved.csv');
    NewFolderPath = strcat('', pwd, '\StatResultsClustering', char(trial), 'time');
    ScriptPath = strcat('', pwd, '\BoxTestAndANOVA_FC.R'); % R script path
    
    % Run R Script of Analysis
    threshold = 'FALSE'; % Needed to run correctly the R code
    command = ['"' 'C:/PROGRA~1/R/R-43~1.0/bin/Rscript" "' ScriptPath ...
        '" "' DataGlobalOutliers '" "' GroupBy '" "' Variable1 '" "' Id '" "' NewFolderPath ...
        '" "' threshold '"'];
    system(command);
    
    movefile("DataGlobalOutliersRemoved.csv", NewFolderPath);
    movefile("OutliarsTable.csv", NewFolderPath);

    %%

    c = 1;
    for gr = unique(Clustering_data.Sex).'
        F_FillArea(means(c, :), ...
            sds(c, :)./sqrt(ns(c)), ...
            cell2mat(Experiment.Project.Palette(gr)), ...
            1:length(means))
        hold on
        c = c+1;
    end
    
    c = 1;
    for gr = unique(Clustering_data.Sex).'
        plot(means(c, :), "LineWidth", 2, 'Color', ...
            cell2mat(Experiment.Project.Palette(gr)))
        hold on
        c = c+1;
    end
    
    ylabel("% in high-activity supercluster")
    xticks([1:length(means)])
    xticklabels(Epochs)
    legend([unique(Clustering_data.Sex).', "", ""])
    hold off
    box off
exportgraphics(gcf, fig_filename, "Append", true, "ContentType", 'vector')    
%% Performing herarchical clustering for the non-normalised data
close all
% Normalising the centroids
% Computing distances between the centroids
Distances_amp = pdist(Centroids); 
% And generating a tree given proximity
Tree = linkage(Distances_amp, "ward"); 

% Generating clustering given the herarchy tree
H_Clusters_amp = cluster(Tree, 'maxclust', 2);
leafOrder_2 = optimalleaforder(Tree, Distances_amp);
cutoff = median([Tree(end-1,3) Tree(end,3)]); % Setting the cutoff

Clustering_data.Supercluster_AMP = ...
    H_Clusters_amp(Clustering_data.Clusters);


% Removing the abnormal cluster
Removed_13_cents = Centroids([1:12, 14:end], :);
Removed_13_cents_Ks = H_Clusters_amp([1:12, 14:end], :);

% For consistency
k_order = unique(H_Clusters_amp(leafOrder_2));

spt_a = F_GetSptIxes(1:5, 31, [1:7]);
subplot(10, 31, spt_a)
    dendrogram(Tree, 'ColorThreshold', cutoff, "Reorder", leafOrder_2);
    hold on
    yline(cutoff, "LineWidth", 2, "Color", 'r', "LineStyle", ":")
    xlabel("Cluster")
    ylabel("Ward distance")


% Setting to equal ylims
yl = [];
spt_b = F_GetSptIxes(7:10, 31, [1:3]);
subplot(10, 31, spt_b)
    area([1:length(Centroids)]./30, ...
        mean(Removed_13_cents(Removed_13_cents_Ks == k_order(1), :)))
    xlabel("Time")
    ylabel("Mean of centroids")
    yl = [yl, ylim()];

spt_c = F_GetSptIxes(7:10, 31, [5:7]);
subplot(10, 31, spt_c)
    xlabel("Time")
    area([1:length(Centroids)]./30, ...
        mean(Removed_13_cents(Removed_13_cents_Ks == k_order(2), :)))
    yl = [yl, ylim()];

subplot(10, 31, spt_b)
    ylim([min(yl, [], "all"), max(yl, [], "all")])
    box off    
subplot(10, 31, spt_c)
    ylim([min(yl, [], "all"), max(yl, [], "all")])
    box off
% Generating indecies
spt_d = F_GetSptIxes(1:10, 31, [9:16])
subplot(10, 31, spt_d)
    Frequencies = F_FrequencyTable(Clustering_data.Epoch,...
           Clustering_data.Supercluster_AMP, true);

    bar(table2array(Frequencies(:, 2:end)).', 'stacked')
        % Labelling it
    xticklabels(Frequencies.Properties.VariableNames(2:end))
    xlim([.5, size(Frequencies, 2)-1+.5])
    ylim([0, 100])
    box off    
    ylabel("% of cases")


spt_e = F_GetSptIxes(1:10, 31, [18:21])
subplot(10, 31, spt_e)
    Frequencies = F_FrequencyTable(Clustering_data.Sex,...
           Clustering_data.Supercluster_AMP, true);
    bar(table2array(Frequencies(:, 2:end)).', 'stacked')
            % Labelling it
    xticklabels(Frequencies.Properties.VariableNames(2:end))
    xlim([.5, size(Frequencies, 2)-1+.5])
    ylim([0, 100])
    box off


%%
spt_f = F_GetSptIxes(1:10, 31, [23:31])
subplot(10, 31, spt_f)

    % Looking at it on an animal-animal basis
    SexEpoch_F = {};
    Frec_ = [];
    Animal = unique(Clustering_data.Animal);
    Sex = [];
    means = [];
    sds = [];
    ns = [];
    c = 1;
    for gr = unique(Clustering_data.Sex).'
        gr_data = Clustering_data(Clustering_data.Sex == gr, :);
        gr_f = [];
        for an = unique(gr_data.Animal).'
            an_data = gr_data(gr_data.Animal == an, :);
            Frequency = table2array(F_FrequencyTable(an_data.Epoch, ...
                an_data.Supercluster_AMP, true))
            gr_f = [gr_f; Frequency(2, 2:end)]
        end
        Sex = [Sex; repelem(gr, length(unique(gr_data.Animal))).'];
        SexEpoch_F{c} = gr_f
        c = c+1;
        means = [means; mean(gr_f, 1)];
        sds = [sds; std(gr_f, [], 1)];
        ns = [ns, size(gr_f, 1)];
        Frec_ = [Frec_; gr_f];
    end
    Frec_;

    TableForStats = table();
    TableForStats.Sex = Sex;
    TableForStats.Animal = Animal;
    TableForStats.Probability = Frec_
    c = 1;

exportgraphics(gcf, fig_filename, "Append", true, "ContentType", 'vector')
    %%

    writetable(TableForStats, "TableForStats.csv");

    % Variables For Grubbs Test
    GroupBy = ['Sex'];
    Variable1 = ['Probability'];
    ScriptPath = strcat('', pwd, '\GrubbsTestOutliersMock.R'); % R script path
    DataTable = char("TableForStats.csv")
    
    % Run R Script of Analysis
    threshold = 'FALSE'; % Needed to run correctly the R code
    command = ['"' 'C:/PROGRA~1/R/R-43~1.0/bin/Rscript" "' ScriptPath ...
        '" "' DataTable '" "' GroupBy '" "' Variable1 '" "' threshold '"'];
    system(command);
    
    % Outliers Table 
    Output = readtable("OutliarsTable.csv");
    
    DataTable2 = readtable(DataTable);
    DataTable2.Animal = string(DataTable2.Animal);
    DataTable2.Sex = string(DataTable2.Sex);
    
    % To Remove Outliers 
    if width(Output) ~= 1;
        for x = height(Output):-1:1;
            outlier = string(table2array(Output(x,2)));
            [index,~] = find(DataTable2{:,:}==outlier);
            DataTable2(index,:) = []; % deletion
        end 
    end
    
    writetable(DataTable2, "DataGlobalOutliersRemoved_2.csv");
    
    % Variables For Statistical Analysis
    Id = ['Animal'];
    DataGlobalOutliers = strcat('', pwd, '\DataGlobalOutliersRemoved_2.csv');
    NewFolderPath = strcat('', pwd, '\StatResultsClustering_', char(trial), 'amp');
    ScriptPath = strcat('', pwd, '\BoxTestAndANOVA_FC.R'); % R script path
    
    % Run R Script of Analysis
    threshold = 'FALSE'; % Needed to run correctly the R code
    command = ['"' 'C:/PROGRA~1/R/R-43~1.0/bin/Rscript" "' ScriptPath ...
        '" "' DataGlobalOutliers '" "' GroupBy '" "' Variable1 '" "' Id '" "' NewFolderPath ...
        '" "' threshold '"'];
    system(command);
    
    movefile("DataGlobalOutliersRemoved_2.csv", NewFolderPath);
    movefile("OutliarsTable.csv", NewFolderPath);
    %%
    for gr = unique(Clustering_data.Sex).'
        F_FillArea(means(c, :), ...
            sds(c, :)./sqrt(ns(c)), ...
            cell2mat(Experiment.Project.Palette(gr)), 1:length(means))
        hold on
        c = c+1;
    end
    
    c = 1;
    for gr = unique(Clustering_data.Sex).'
        plot(means(c, :), "LineWidth", 2, 'Color', ...
            cell2mat(Experiment.Project.Palette(gr)))
        hold on
        c = c+1;
    end
    
    labels_ = Epochs;
    ylabel("% in high-activity supercluster")
    xticks([1:length(means)])
    xticklabels(labels_)
    legend([unique(Clustering_data.Sex).', "", ""])
    hold off
    box off

f_ = gcf;
f_.Position = [312,449,1600,530];
exportgraphics(f_, "AmplitudeSuperclustering_GLOBAL.pdf", ...
    'ContentType', 'vector')
exportgraphics(gcf, fig_filename, "Append", true, "ContentType", 'vector')

%%
y = reshape(Frec_, 1, []).';
x = reshape(repmat(1:15, 13, 1), 1, []).';
t = table();
t.x = x; t.y = y;


mdl = fitlm(t);
plot(mdl)
hold on
scatter(x, y, 27, 'k', 'filled');

[R,P] = corrcoef(x, y)
%% GENERATING TABLE FOR THE LGMM