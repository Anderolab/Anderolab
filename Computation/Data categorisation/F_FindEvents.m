function[Events_Merg, Events] = F_FindEvents(Data, LowBound, ...
    HighBound, Visualise)
%% Determining the frames of interest
% Max within high activity frames
max_ = find(islocalmax(Data(1, :)).*(Data(1, :) > HighBound) == 1);
% Min within low activity frames
min_ = find(islocalmin(Data(1, :)).*(Data(1, :) < LowBound) == 1);

%% Identifying high activity epochs
% For indexing
start_ = 1;

% Storage array
Events = [];

% For visualisation purposes be that needed
area_ = zeros(1, length(Data));

% Computing the epoch start and end frames
while true
    % Finding the next min
    end_ = find(min_ > (max_(start_)), 1);

    % For visualisation purposes be that needed
    area_(max_(start_):min_(end_)) = 1;

    % Ending the loop at the end if it finishes in a high activity period
    if isempty(end_)
        break
    end

    % Saving the event
    Events = [Events; max_(start_), min_(end_)];
    
    % Identifying the next start of epoch
    start_ = find(max_ > (min_(end_)), 1);

    % Ending the loop at the end if it finishes in a low activity period
    if isempty(start_)
        break
    end
end
Events
plot(Data)
hold on
yline(LowBound)
yline(HighBound)

%% Merging events that are close by
% Storage array
Events_Merg = Events(1, :);

% For indexing purpouses
c = 2;

% Looping through each event
for event = 2:size(Events, 1)

    % Concatenating if they're not distant
    if Events_Merg(c-1, 2) >= (Events(event, 1)-10)
        Events_Merg(c-1, 2) = Events(event, 2);
    else

        % Saving when they are
        Events_Merg(c, :) = Events(event, :);
        c = c+1;
    end

end
if Visualise
    merg_area = zeros(1, size(Data, 2));
    for i = 1:size(Events_Merg, 1)
        merg_area(Events_Merg(i, 1):Events_Merg(i, 2)) = 1;
    end

    area((Data+abs(min(Data))).*merg_area, 'FaceColor', 'k', ...
        'FaceAlpha', 1)
    hold on
    yline(LowBound+abs(min(Data)), 'Color', 'r', 'LineWidth', 2)
    yline(HighBound+abs(min(Data)), 'Color', 'r', 'LineWidth', 2)
    plot(Data+abs(min(Data)), 'Color', 'k')
    yticks([])
    ylabel('Movement score (u)')
    xlabel('Time (frames)')

end
end