function[] = F_HeatScatter(x, y, bins)
%F_HEATSCATTER Summary of this function goes here
%   Detailed explanation goes here

% Setting the limits for the distribution
min_x = min(x);
min_y = min(y);
max_x = max(x);
max_y = max(y);

% Setting an edge to the data
edge_x = abs(max_x - min_x)*.05;
edge_y = abs(max_y - min_y)*.05;

% Creating bin limits
pts_x = linspace(min_x-edge_x, max_x + edge_x, bins);
pts_y = linspace(min_y-edge_y, max_y + edge_y, bins);

% Calculating density per bin
N = histcounts2(y, x, pts_y, pts_x);



% Plotting
imagesc(pts_x, pts_y, N);
hold on
set(gca, 'XLim', pts_x([1, end]), 'YLim', pts_y([1, end]), 'YDir', 'normal');
hold off
end

