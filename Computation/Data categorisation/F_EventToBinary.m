function [Boolean] = F_EventToBinary(Dataset, Event)
    Boolean = zeros(1, size(Dataset, 2));

    % Looping per event
    for event = 1:size(Event, 1)
        Boolean(1, Event(event, 1):Event(event, 2)) = 1;
    end
    area(Boolean.*Dataset)
    hold on
    plot(Dataset)
    hold off
end

