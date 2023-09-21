
function [Peaks] = F_ExtractFrequency(dataset, binwidth,recording_frequency)

peak_results = zeros(size(dataset));
raw_peak_results = peak_results;
frequency_results = cat(2, size(dataset), 1);
frequency_results(2) = floor(frequency_results(2)/binwidth);
frequency_results = zeros(frequency_results);

% Looping throug the trials and neurons
for trial = 1:size(dataset, 3)
    for neuron = 1:size(dataset, 1)

        % Extracting the peaks %% Vaidate the interval
        [raw_peaks, peaks] = F_FilterPeaks(dataset(neuron, :, trial), 40);

        % Saving
        peak_results(neuron, :, trial) = peaks;
        raw_peak_results(neuron, :, trial) = raw_peaks;


        % Computing the frequency
        for bin = 1:size(frequency_results, 2)
            spikes_in_bin = sum(peaks(:, 1+(bin-1)*binwidth:(bin*binwidth)));
            frequency = spikes_in_bin/(binwidth/recording_frequency);
            frequency_results(neuron, bin, trial) = frequency;
        end
    end
end

Peaks.frequency = frequency_results;
Peaks.peaks = peak_results;
Peaks.binned = double(frequency_results > 0);
Peaks.raw_peaks = raw_peak_results;

end