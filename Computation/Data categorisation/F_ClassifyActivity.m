%% Getting the data
pks1 = IMO20.Cat.Dat_PEAKS;
pks2 = IMO10.Cat.Dat_PEAKS;

%% Computing frequency
% Unified array.
Data = {pks1, pks2};

% Output storage and counter
Frequency = {};
fs = [];
c_ = 1;

% Frequency calculation
for epoch = Data
    Frequency{c_} = sum(epoch{:}, 2)/(size(epoch{:}, 2)/RF);
    fs = [fs; Frequency{c_}];
    c_ = c_ + 1;
end

%% Computing AUC
flu1 = IMO20.Cat.Dat_ALLN;
flu2 = IMO10.Cat.Dat_ALLN;

Data = {flu1, flu2};

% Output storage and counter
Fluorescence = {};
flu = [];
c_ = 1;

% Frequency calculation
for epoch = Data
    Fluorescence{c_} = sum(epoch{:}, 2)/(size(epoch{:}, 2)/RF);
    flu = [flu; Fluorescence{c_}];
    c_ = c_ + 1;
end
%% Creating a classifier

fit_data = fs;
% Generating histogram
hg = histogram(fit_data , 50, "EdgeColor", "k", "FaceColor", "w",...
    "Normalization", 'probability');
hold on

% Getting centre of bins and fraction of total
bin_cent = hg.BinEdges(1:end-1) + hg.BinWidth(1)/2;

% Performing the fit
f = fit(bin_cent.', hg.Values.',  'gauss2');

% Getting centre of bins and fraction of total
bin_cent = hg.BinEdges(1:end-1) + hg.BinWidth(1)/2;

% Visualising individual points
scatter(bin_cent, hg.Values, 5, 'k', 'filled')
hold on

area_ = f.b1 + 1.96*(f.c1/sqrt(2)):.002:max(bin_cent);
f_s = [];
for i = area_
    f_s = [f_s, f(i)];
end
area(area_, f_s, 'FaceAlpha', .3)



% Plotting the fit
f_plot = plot(f);


% Visualisation
set(f_plot, "LineWidth", 2)
legend(["", "", "Tonic", "1g fit"])
ylabel('Probability')
xlabel('Spiking requency')
xlim([0, max(bin_cent)+ hg.BinWidth(1)/2])
hold off

% xline(f.b1)
% xline(f.b1 + 1.96*(f.c1/sqrt(2)))
% hold off