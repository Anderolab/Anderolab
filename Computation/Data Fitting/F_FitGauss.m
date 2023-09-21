function [f] = F_FitGauss(fit_data, model)

% Generating histogram
hg = histogram(fit_data , "BinMethod", 'fd', "EdgeColor", "k", "FaceColor", "w",...
    "Normalization", 'probability');
hold on
reshaped = reshape(fit_data, 1, []);
xl = xlim();
yl = ylim();
%
% Getting centre of bins and fraction of total
bin_cent = hg.BinEdges(1:end-1) + hg.BinWidth(1)/2;

% Performing the fit
f = fit(bin_cent.', hg.Values.',  model);

% Getting centre of bins and fraction of total
bin_cent = hg.BinEdges(1:end-1) + hg.BinWidth(1)/2;

% Visualising individual points
scatter(bin_cent, hg.Values, 5, 'k', 'filled')
hold on

% For clarity
legnd = ["", ""];

if sum(model ~= 'gauss1')

    % Plotting individual gauss
    coffs = coeffvalues(f);
    for coff_ix = 1:length(coffs)/3

        % Plotting
        plot(-10:.001:1, F_ModGauss(-10:.001:1, coffs((coff_ix-1)*3+1), ...
            coffs((coff_ix-1)*3+2), coffs((coff_ix-1)*3+3)), ...
            "LineWidth", 3);
        hold on

        legnd = [legnd, strcat("Gauss", num2str(coff_ix))];
    end
    %area(area_, f_s, 'FaceAlpha', .3)
end


% Plotting the fit
f_plot = plot(f);
legnd = [legnd, "Fit"];


ylim(yl)
% Visualisation
set(f_plot, "LineWidth", 2)
legend(legnd)
ylabel('Probability')
xlim(nanmean(reshaped) + 3*[-nanstd(reshaped), nanstd(reshaped)])
hold off
end

%% 
