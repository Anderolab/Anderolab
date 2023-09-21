%% Loading pre-requisites
IMOE = load('All_IMO_Data_E.mat');
IMOL = load('All_IMO_Data_L.mat');


%% Saving structures
StruggleScore = []; % SDs away from boostrapping noise
AIX = []; % Animal number corresponding to each neuron
NIX = []; % Neuron index within the animal
Means = [];
SDs = [];
% Iterating each animal
for animal = IMOE.IMO_Data.Animals

    % To binary
    Event = F_EventToBinary(Data(IMOE.IMO_Data.Animals == animal, :), ...
        Events{1}{IMOE.IMO_Data.Animals == animal});
    
    % Selecting the animal relevant dataset
    selected_data = ...
        IMOE.IMO_Data.RawTraces(IMOE.IMO_Data.AIX == animal, :);
    selected_data = selected_data(:, ~isnan(selected_data(1, :)));
    selected_data = ... % DeltaF given trace median
        (selected_data-prctile(selected_data, 5, 2)) ./ ...
        prctile(selected_data, 5, 2); % Removal of empty frames
    


    Event = Event(:, ~isnan(selected_data(1, :))); % Empty frame removal
    

    % Computing the score
    R_L = nanmean(selected_data(:, Event == 1), 2); % During movement
    R_S = nanmean(selected_data(:, Event == 0), 2); % During stationary
    Score = (R_L-R_S)./(R_L+R_S); % Computing the score


    % Generating the boostrap save
    Bootstrap_results = zeros(size(selected_data, 1), 1000);
    
    % Bootstrapping
    for i = 1:1000

        % Generaring randomised index
        ix = randperm(size(Event, 2));
        rand_data = selected_data(:, ix); % Randomising

        % Computing score for randomised dataset
        R_L = nanmean(rand_data(:, Event == 1), 2); % Movement
        R_S = nanmean(rand_data(:, Event == 0), 2); % Stationery
        Bootstrap_results(:, i) = (R_L-R_S)./(R_L+R_S); % Saving scores
    end

    
    % Expressing struggle score as SDs away from mean
    SDs = [SDs; nanstd(Bootstrap_results, [], 2)];
    Means = [Means; mean(Bootstrap_results, 2)];
    StruggleScore = [StruggleScore; Score];
    % Generating neuron index and animal index saving array
    NIX = [NIX, 1:size(selected_data, 1)];
    AIX = [AIX, repelem(animal, size(selected_data, 1))];
end

%% Visualisation

% Creating an ideal event
Pre_Frames = 40;
Post_Frames = 70;

% Generating the pattern
pattern = [repelem(0, Pre_Frames), repelem(1, Post_Frames)];

% Identifying the animal and neuron to visualise
animal = AIX(find(StruggleScore  == max(StruggleScore), 1));
n_ix = NIX(find(StruggleScore == max(StruggleScore), 1));

% Extracting the animal-relevant event binary
Event = F_EventToBinary(Data1(IMOE.IMO_Data.Animals == animal, :), ...
        Events{1}{IMOE.IMO_Data.Animals == animal});

% Extracting the Calcium trace for the relevant neuron
selected_data = IMOE.IMO_Data.RawTraces(IMOE.IMO_Data.AIX == animal, :);
selected_data = (selected_data-nanmedian(selected_data, 2)) ./ ...
    nanmedian(selected_data, 2);

% Identifying timepoints that meet the ideal event criterion
finds = strfind(Event, pattern);

close all

% Generating storage arrays for the visualisation
store_mov = []; % Movement
nn = []; % Neuronal signal

for find_ = finds
    % Getting pattern-associated movement trace
    store_mov = [store_mov; Data1(IMOE.IMO_Data.Animals == animal, ...
        find_:(find_+Pre_Frames + Post_Frames))];
    % Getting the calcium trace associated to the ideal event
    nn = [nn; selected_data(n_ix, find_:(find_+Pre_Frames + Post_Frames))];
end

% Actual visualisation
subplot(2, 2, 1)
xline(Pre_Frames)
xlabel('Time (frames)')
ylabel('Struggle score')
hold on
F_FillArea(mean(store_mov, 1), ...
    std(store_mov, 1)./sqrt(size(store_mov, 1)), [0, 75, 75]./255, ...
    1:size(store_mov, 2))
plot(mean(store_mov, 1), 'Color', [0, 75, 75]./255)
hold off

subplot(2, 2, 3)
F_FillArea(wdenoise(mean(nn, 1)), std(nn, 1)./sqrt(size(nn, 1)), ...
    [0, 75, 75]./255, 1:length(mean(nn, 1)))
hold on
xline(Pre_Frames)
plot(wdenoise(mean(nn, 1)), 'Color', [0, 75, 75]./255)
hold off
xlabel('Time (frames)')
ylabel('\DeltaF/F_0')

s_ = subplot(2, 2, 2);
plot(Data1(IMOE.IMO_Data.Animals == animal, :))

s_s = subplot(2, 2, 4);
area(selected_data(n_ix, :))
linkaxes([s_, s_s], 'x')