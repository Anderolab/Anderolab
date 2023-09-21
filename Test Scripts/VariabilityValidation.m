% Required variables
Variable = Centroids;
grs = H_Clusters_amp;
ReportName = 'Centroids';
Time = 1:length(Variable);
savename = [];

ReportName = strcat(ReportName, '.pdf');

%% STEP 1 - BINNING SINCE NEED FOR LESS DIMENSIONS
Bins = size(Variable, 1)-1;
BW  = discretize(Time, Bins);
Binned = zeros(size(Variable, 1), Bins);

for bin = 1:Bins
    Binned(:, bin) = mean(Variable(:, BW == bin), 2);
end

% Viewing the centroids
for centroid = 1:size(Variable, 1)
    plot(Binned(centroid, :) + centroid);
    hold on
end
hold off
close all
%% STEP 2 - Performing a PCA
[~, score, ~, ~, explained] = pca(Binned);
CF = cumsum(explained);
MaxPCs = find(CF>99.9, 1);
MinPCs = find(CF>95, 1);
xline(MinPCs, "LineWidth", 2, "Color", 'r', 'LineStyle', ":")
hold on
yline(CF(MinPCs), "LineWidth", 2, "Color", 'r', 'LineStyle', ":")
plot(CF(1:MaxPCs), "Color", 'k', 'LineWidth', 2);
bar(explained(1:MaxPCs), 'FaceColor', 'w')
bar(explained(1:MinPCs), 'FaceColor', 'r')
hold off
ylabel('% of variability explained')
xlabel("PC")
legend(["", "95% Threshold", 'Cumulative', 'Single PC'])
f = gcf;
f.Position = [680,558,382,420];
exportgraphics(f, ReportName, "Append", true, "ContentType", 'Vector')
% close all
%% STEP 3 - Gathering the variables
AUC = trapz(Binned, 2); % Area under the curve
TOM = []; % Time of maximum fluorescence
MaxF = max(Binned, [], 2);
MinF = min(Binned, [], 2);

% Time of max fluorescence
for centroid = 1:1:size(Variable, 1)
    TOM = [TOM; find(Binned(centroid, :) == MaxF(centroid))];
end

% Setting the variable names and variable array
VarNames = ["AUC", "TOM", "MaxF", "MinF"];
Vars = table(AUC, TOM, MaxF, MinF, 'VariableNames', VarNames);

PCNames = "PC" + string(1:MinPCs);
PCs = table();
c = 1;
for pc = PCNames
    PCs.(pc) = score(:, c);
    c = c+1;
end
%% Performing the analysis

c = 1;
R_Scores = zeros(length(VarNames), length(PCNames));
Acc_Scores = R_Scores;
for i = 1:length(VarNames)
    for j = 1:length(PCNames)

        % Setting the subplot
        subplot(length(VarNames), length(PCNames), c)

        % Generating a linear regression model
        mdl = fitlm(PCs.(PCNames(j)), Vars.(VarNames(i)));
        plot(mdl)
        hold on
        gscatter(PCs.(PCNames(j)), Vars.(VarNames(i)), grs)

        % For consistency
        g = gcf;
        xlabel("")
        ylabel("")
        xl = xlim();
        yl = ylim();

        if j == 1
            ylabel(VarNames(i))
        end

        if i == length(VarNames)
            xlabel(strcat(PCNames(j), "; ", ...
                num2str(round(explained(j), 3)),' var'))
        end
        title("")

        if mdl.Coefficients.pValue(2) <= .05
            set(gca,'Color', [255, 215, 212]./255)

            title(strcat("p = ", ...
            num2str(round(mdl.Coefficients.pValue(2), 5))), ...
            strcat("R^2: ",...
            num2str(round(mdl.Rsquared.Adjusted, 3))));
        end
        
        % Generating an LDA for validation
        MdlLinear = ...
            fitcdiscr([PCs.(PCNames(j)), Vars.(VarNames(i))], grs);
        K = MdlLinear.Coeffs(1, 2).Const;  
        L = MdlLinear.Coeffs(1, 2).Linear;
        f = @(x1,x2) K + L(1)*x1 + L(2)*x2;

        % Visualising the boundary
        hold on
        h2 = fimplicit(f,[min(xl), max(xl), min(yl), max(yl)]);
        xlim(xl)
        ylim(yl)
        h2.Color = 'k';
        h2.LineWidth = 2;
        
        % Testing the accuracy of the LDA
        test_ = ...
            predict(MdlLinear,[PCs.(PCNames(j)), Vars.(VarNames(i))]);
        cf = confusionmat(grs, test_);
        accuracy = sum(diag(cf), 'all')/sum(cf, 'all');
        
        % Saving
        R_Scores(i, j) = mdl.Rsquared.Adjusted;
        Acc_Scores(i, j) = accuracy;
        c = c+1;
        l = legend([""]);
        set(l,'visible','off')

    end
end
sgtitle(strcat("df = ", num2str(mdl.DFE)))
f = gcf;
f.Position = [680,209,917,769];
% exportgraphics(f, ReportName, "Append", true, "ContentType", 'Vector')
% close all
%%
% Last visualisation
imagesc(R_Scores)
yticks(1:length(VarNames))
xticks(1:length(PCNames))
colorbar
xticklabels(PCNames)
yticklabels(VarNames)
title("R scores")
f = gcf;
f.Position = [680,209,500,420];
exportgraphics(f, ReportName, "Append", true, "ContentType", 'Vector')
close all

%% Visualising accuracy
imagesc(Acc_Scores)
yticks(1:length(VarNames))
xticks(1:length(PCNames))
colorbar
xticklabels(PCNames)
yticklabels(VarNames)
title("R scores")
f = gcf;
f.Position = [680,209,500,420];
exportgraphics(f, ReportName, "Append", true, "ContentType", 'Vector')
close all

%% Testing reliability of the model
% Generating synthetic y
rng(1)
y = randn(1, size(Variable, 1)).';

observed = [];
chance = [];
sd = [];
test_diff = [];
test_stats = {};

spt_cols = length(PCNames);
for pc = 1:length(PCNames)
    x = score(:, pc);
     % Generating an LDA for validation
    MdlLinear = ...
        fitcdiscr([x, y], grs);
    K = MdlLinear.Coeffs(1, 2).Const;  
    L = MdlLinear.Coeffs(1, 2).Linear;
    f = @(x1,x2) K + L(1)*x1 + L(2)*x2;
    
    % Testing the accuracy of the LDA
    test_ = ...
        predict(MdlLinear,[x, y]);
    cf = confusionmat(grs, test_);
    subplot(5, spt_cols, pc+4*spt_cols)
    imagesc(cf)
    accuracy = sum(diag(cf), 'all')/sum(cf, 'all');
    observed = [observed, accuracy];
    subplot(5, spt_cols, pc+3*spt_cols)
    gscatter(x, y, grs)
    hold on
    xl = xlim();
    yl = ylim();

    l = legend([""]);
    set(l,'visible','off')
    h2 = fimplicit(f,[min(xl), max(xl), min(yl), max(yl)]);
    xlim(xl)
    ylim(yl)
    h2.Color = 'k';
    h2.LineWidth = 2;
    hold off
    shuffled = [];
    % Now onto the randomising
    for i = 1:1000
        rand_grs = grs(randperm(length(grs)));

        MdlLinear = ...
        fitcdiscr([x, y], rand_grs);
        K = MdlLinear.Coeffs(1, 2).Const;  
        L = MdlLinear.Coeffs(1, 2).Linear;
        f = @(x1,x2) K + L(1)*x1 + L(2)*x2;
        
        % Testing the accuracy of the LDA
        test_ = ...
            predict(MdlLinear,[x, y]);
        cf_ = confusionmat(rand_grs, test_);
        shuffled = [shuffled, sum(diag(cf_), 'all')/sum(cf_, 'all')];

    end
    chance = [chance, mean(shuffled)];
    [h, p, ~, stat] = ttest(shuffled, accuracy);
    test_diff = [test_diff; [h, p]];
    test_stats{pc} = stat;
    sd = [sd, std(shuffled)];


end
subplot(5, spt_cols, 1:(spt_cols*3))
F_FillArea(chance, sd*1.96, 'k', 1:length(PCNames));
xlim([.5, length(PCNames)+.5])
xticks(1:length(PCNames))
xticklabels(PCNames)
ylabel("Model performance")
hold on
plot(chance, 'Color', 'k', 'LineWidth', 2)
plot(observed, 'Color', 'r', 'LineWidth', 2)

Z = (observed-chance)./sd;

%
format longg
P = @(Z) erfc(-Z/sqrt(2))
P(-abs(Z))

%% Testing reliability of the model
% Generating synthetic y
rng(1)
y = randn(1, size(Variable, 1)).';

observed = [];
chance = [];
sd = [];
test_diff = [];
test_stats = {};

spt_cols = length(VarNames);
for pc = 1:length(VarNames)
    x = Vars.(VarNames(pc));
     % Generating an LDA for validation
    MdlLinear = ...
        fitcdiscr([x, y], grs);
    K = MdlLinear.Coeffs(1, 2).Const;  
    L = MdlLinear.Coeffs(1, 2).Linear;
    f = @(x1,x2) K + L(1)*x1 + L(2)*x2;
    
    % Testing the accuracy of the LDA
    test_ = ...
        predict(MdlLinear,[x, y]);
    cf = confusionmat(grs, test_);
    subplot(5, spt_cols, pc+4*spt_cols)
    imagesc(cf)
    accuracy = sum(diag(cf), 'all')/sum(cf, 'all');
    observed = [observed, accuracy];
    subplot(5, spt_cols, pc+3*spt_cols)
    gscatter(x, y, grs)
    hold on
    xl = xlim();
    yl = ylim();

    l = legend([""]);
    set(l,'visible','off')
    h2 = fimplicit(f,[min(xl), max(xl), min(yl), max(yl)]);
    xlim(xl)
    ylim(yl)
    h2.Color = 'k';
    h2.LineWidth = 2;
    hold off
    shuffled = [];
    % Now onto the randomising
    for i = 1:1000
        rand_grs = grs(randperm(length(grs)));

        MdlLinear = ...
        fitcdiscr([x, y], rand_grs);
        K = MdlLinear.Coeffs(1, 2).Const;  
        L = MdlLinear.Coeffs(1, 2).Linear;
        f = @(x1,x2) K + L(1)*x1 + L(2)*x2;
        
        % Testing the accuracy of the LDA
        test_ = ...
            predict(MdlLinear,[x, y]);
        cf_ = confusionmat(rand_grs, test_);
        shuffled = [shuffled, sum(diag(cf_), 'all')/sum(cf_, 'all')];

    end
    chance = [chance, mean(shuffled)];
    [h, p, ~, stat] = ttest(shuffled, accuracy);
    test_diff = [test_diff; [h, p]];
    test_stats{pc} = stat;
    sd = [sd, std(shuffled)];


end
subplot(5, spt_cols, 1:(spt_cols*3))
F_FillArea(chance, sd*1.96, 'k', 1:length(VarNames));
xlim([.5, length(VarNames)+.5])
xticks(1:length(VarNames))
xticklabels(VarNames)
ylabel("Model performance")
hold on
plot(chance, 'Color', 'k', 'LineWidth', 2)
plot(observed, 'Color', 'r', 'LineWidth', 2)

Z = (observed-chance)./sd;

%
format longg
P = @(Z) erfc(-Z/sqrt(2))
P(-abs(Z))
