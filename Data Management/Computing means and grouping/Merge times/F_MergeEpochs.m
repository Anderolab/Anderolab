function [Output] = F_MergeEpochs(SourceTable, Variable, Episodes,...
    Cathegory, Task, Pre_Frames, Post_Frames)
% Temporary storage of the tables
tables = {};
Cat = [];

% Looping through the columns
c = 1;  % Counter function
for ser = Episodes.'
    tables{c} = F_SelectEpoch(SourceTable, Variable, ser.', Task, ...
        Pre_Frames, Post_Frames); % Splits
    Cat = [Cat, repelem(Cathegory(c), size(tables{c}, 1))];
    c = c+1;
end
Output = vertcat(tables{:}); % Merges
Output.Cathegory = Cat.';

% Sorting
[~, idx] = ismember(Output.Cathegory, unique(Cathegory, 'stable'));
[~, sortorder] = sort(idx);
Output = Output(sortorder,:);

end

