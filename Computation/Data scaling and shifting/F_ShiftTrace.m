function [Output] = F_ShiftTrace(SourceTable, Variable, Epoch, Method, ...
    Task, Divide)
% Gathering the source to use as the shift

if isempty(Epoch) == false

    Ref = F_SelectEpoch(SourceTable, Variable, Epoch, Task, 0, 0);

else
    "H"
    Ref = SourceTable;
end

if Method == "Mean"
    shift_ = mean(Ref.(Variable), 2)
elseif contains(Method, "Pctile")
    ptile = split(Method, " ");
    shift_ = prctile(Ref.(Variable), double(ptile(2)), 2);
    clear ptile
end
% Computing the mean

if Divide
    SourceTable.(Variable) = (SourceTable.(Variable) - shift_)./abs(shift_);
else
    SourceTable.(Variable) = SourceTable.(Variable) - shift_;
end

% Saving the output
Output = SourceTable;
Output.Shift = -mean(Ref.(Variable), 2);
end

