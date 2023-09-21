function [Concatenated] = F_ConcatenateField(Data, Field)
%F_CONCATENATEFIELD Summary of this function goes here
%   Detailed explanation goes here
an_s = fieldnames(Data);
Animals = string(regexp(string(an_s(1:end-1)), '\d*', 'Match'));

% Pre-computing length of save struct
counter_neurons = zeros(size(Animals));
counter_frames = zeros(size(Animals));
c = 1;
for animal = Animals.'
    counter_neurons(c) = size(Data.("M" + animal).(Field), 1);
    counter_frames(c) = size(Data.("M" + animal).(Field), 2);
    c = c+1;
end

% Setting an index
ix_ = [cumsum([1, counter_neurons(1:end-1).']); ...
    cumsum([1, counter_neurons(1:end-1).']) + counter_neurons.'-1].';
% And max number of neurons
neurons = sum(counter_neurons, 'all');
% Generating the empty save table
Concatenated = table();
Concatenated.Animal = repelem("", neurons).';
Concatenated.Neuron = zeros(neurons, 1);
Concatenated.(Field) = zeros(neurons, min(counter_frames));

% Populating it
c = 1;
for animal = Animals.'
    Concatenated.(Field)(ix_(c, 1):ix_(c, 2), :) = ...
        Data.("M" + animal).(Field)(:, 1:min(counter_frames));
    Concatenated.Animal(ix_(c, 1):ix_(c, 2), :) = ...
        repelem(animal, counter_neurons(c)).';
    Concatenated.Neuron(ix_(c, 1):ix_(c, 2), :) = ...
        [1:counter_neurons(c)].';
    c = c+1;
end
end

