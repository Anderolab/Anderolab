function [Clustering] = F_PreClustering(Experiment, Trial, Epochs, PreFiltering, p_frames)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%% Alligning all the sessions
an_s = fieldnames(Experiment.(Trial));
Animals = string(regexp(string(an_s(1:end-1)), '\d*', 'Match'));
Report = [];

% Finding all trial lengths
lengths = [];
c = 1;
for animal = Animals.'
    lengths(c) = size(Experiment.(Trial).('M'+animal).("Raw"), 2);
    c = c+1;
end

% Identifying outliers
outliers = isoutlier(lengths, "grubbs");

for i = 1:sum(outliers)
    out_an = Animals(outliers);
    out_len = string(lengths(outliers));
    prompt = strcat("    Animal ", ...
        out_an(i), " was identified as an outlier with ",  ...
        out_len(i), ' frames.');
    fprintf('%s\n', prompt)
    Report = [Report; prompt];
end

boxplot(lengths)

% Croppig all sessions
len = min(lengths);

prompt = strcat("    All sessions will be cropped to ", ...
    num2str(len), " frames.");

fprintf('%s\n', prompt)

Animals = Animals(outliers == 0);

%% Extracting the relevant task
Task = Experiment.(Trial).Task;

%% Generating the save arrays
Clustering = table();
Flu = [];
Animal = [];
Epoch = [];
Neuron = [];
AllNeurons = [];
AllNeurons_r = [];
%% Filtering neuronal responses
if PreFiltering == true

    % GENERATING GLOBALISED ARRAY FOR OUTLIER SELECTION
    for an_ = Animals.'
        AllNeurons = [AllNeurons; Experiment.(Trial).("M"+an_).Filt(:, 1:len)];
        AllNeurons_r = [AllNeurons_r; Experiment.(Trial).("M"+an_).Raw(:, 1:len)];
    end

    % OUTLIER IDENTIFICATION GIVEN THE MEAN FLUORESCENCE
    Outs_MEAN = ...
        F_Outliers(sum(AllNeurons, 2), "median", "Mean flu signal");

    % OUTLIER IDENTIFICATION GIVEN THE PEAKS
    peaks = sum(islocalmax(AllNeurons, 2), 2);
    close all
    plot(peaks)
    Outs_PEAKS = ...
        F_Outliers(sum(peaks, 2), "median", "Number of calcium events");

    % OUTLIER IDENTIFICATION GIVEN THE SUM OF SQUARE DISTANCES
    SSDs = sum((AllNeurons - AllNeurons_r).^2, 2);
    Outs_SSDs = F_Outliers(SSDs, "median", "Sum of square distances");
    size(Outs_SSDs)
    size(Outs_PEAKS)
    size(Outs_MEAN)
    Outs = sum((Outs_SSDs + Outs_PEAKS + Outs_MEAN) > 1);

    close all
    pie([Outs, length(Outs_MEAN)-Outs]);
    legend(["Outliers", ""])


end
%% Gathering all neurons
for an_ = Animals.'
    c_ = 0;
    for epoch = Epochs

        % Setting the start and end
        start = Task.Start(string(Task.Titles) == epoch)-1;
        end_ = Task.Frames(string(Task.Titles) == epoch)+start+p_frames;
        % Selecting the time of interest
        F = Experiment.(Trial).("M"+an_).Filt;
        Flu = [Flu; ...
            (F(:, start:end_)- min(F, [], 2))./(max(F, [], 2)-min(F, [], 2))];
        % Saving the animal, neuron and epoch
        Animal = [Animal; ...
            repelem(an_, size(Experiment.(Trial).("M"+an_).Filt, 1)).'];
        Neuron = [Neuron; [1:size(Experiment.(Trial).("M"+an_).Filt, 1)].']; %#ok<NBRAK1> 
        Epoch = [Epoch; ...
            repelem(epoch, size(Experiment.(Trial).("M"+an_).Filt, 1)).'];
        c_ = c_ +1;
    end
end
Outs = repelem((Outs_SSDs + Outs_PEAKS + Outs_MEAN), c_);
%% Saving the output
Clustering.Flu = Flu(Outs == 0, :);
Clustering.Animal = Animal(Outs == 0);
Clustering.Epoch = Epoch(Outs == 0);
Clustering.Neuron = Neuron(Outs == 0);

end

