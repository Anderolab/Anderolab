function [BestK, BestSSD] = F_GetBestK(Dataset, MaxK)
%F_GETBESTK Summary of this function goes here
%   Detailed explanation goes here

%% STEP 0 - Generating the save array
kscores = [];
Dataset = gpuArray(Dataset); % For speed

%% STEP 1 - Iterating through all the Ks
ktestrange = 1:ceil(MaxK/100):MaxK; % Avoiding +100 iterations
wb = waitbar(0, "Running KMEANs to identify ideal K"); % Notify user
for i = ktestrange % Iterating
    waitbar(i/MaxK) % Notification
    [~, ~, score] = kmeans(Dataset, i, 'MaxIter', 500);

    kscores = [kscores, sum(score)];
    
end
close(wb)
close all

%% STEP 2 - Visualising and computing the knee

subplot(1, 2, 2)
cplot_ = cdfplot(kscores); % cumulative probability transform

x_ = cplot_.XData; % Setting for the visualisation
y_ = cplot_.YData; % ""
plot(x_, y_, "LineWidth", 1.2, "Color", 'k') % Viewing the cumulative prob

box off % Vis
grid off % Vis
hold on % Vis

diag_ = [x_(find(y_ == 0, 1, 'last')), 0; ...
    x_(find(y_ == 1, 1)), 1]; % Points to generate the diagonal function
plot(diag_(:, 1), diag_(:, 2), "LineWidth", 1.2, "Color", 'k', ...
    "LineStyle", "--") % Actual diagonal visualisation

% Transforming the diagonal into a function y = mx+n
    m = (diag_(2, 2)-diag_(1, 2))/(diag_(2, 1)-diag_(1, 1));
    n = diag_(2, 2)-diag_(2, 1)*m;
    y = x_.*m + n;

dists = (y_ - y); % Computing the distances
plot(x_, dists, "LineWidth", 1.2, "Color", 'r'); % And plotting them


BestSSD = x_(dists == max(dists(~isinf(dists)))); % Identify ideal SSDs
% And visualising it
    yline(max(dists(~isinf(dists))), "LineStyle", ":", "LineWidth", 1.2)
    xline(BestSSD, "LineStyle", ":", ...
    "LineWidth", 1.2);

xl = xlim(); % Setting the xlim for consistency with other subplot

xlabel("SSDs") % Vis
ylabel("CDF(SSDs)") % Vis
legend(["CDF", "Diagonal", "Distance to diagonal", "Highest curvature"])
title("Knee identification")
hold off


%% STEP 3 - Visualising the knee location and identifying ideal K
subplot(1, 2, 1)

% Plotting the SSDs given the selected k range
plot(ktestrange, kscores, "LineWidth", 1.2, "Color", 'k')
hold on % Vis

% Computing the ideal k
BestK = ktestrange(kscores == BestSSD);
yline(BestSSD, "LineStyle", ":", "LineWidth", 1.2); % View ideal SSDs
xline(BestK, "LineWidth", 1.2, "Color", 'r'); % View ideal K
scatter(ktestrange(kscores == BestSSD), BestSSD, 30, 'red', ...
    'filled', "LineWidth", 2, "MarkerFaceColor", 'w', ...
    'MarkerEdgeColor', 'r') % Showing knee

xl_ = xlim(); % For text location
fig_txt = ("Knee @ K = " + num2str(ktestrange(kscores == BestSSD)));
text(ktestrange(kscores == BestSSD) + ((xl_(2)-xl_(1))/20), ...
    BestSSD+((xl(2)-xl(1))/20), fig_txt) % Actual labell

legend(["SSDs", "", "", "Knee"]) % Vis
xlabel("K") % Vis
ylabel("SSD(k)") % Vis
box off % Vis
ylim(xl) % Setting limit for consistency purpuses
title("SSD @ different K_s")
hold off

end

