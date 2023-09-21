function [Output] = F_SelectEpoch(SourceTable, Variable, Epoch, ...
    Task, Pre_Frames, Post_Frames)
%F_SELECTEPOCH Summary of this function goes here
%   Detailed explanation goes here

% Getting the variable
Data = SourceTable.(Variable);

% Storage array
Segment = [];

% If buffer frames are wanted
if isempty(Pre_Frames) == false
    ep_ix = string(Task.Titles) == Epoch(1);
    sec = Data(:, ...
        (Task.Start(ep_ix)-Pre_Frames):(Task.Start(ep_ix)-1));
    Segment = [Segment, sec];
end

% Looping per session
for ep = Epoch   
    ep_ix = string(Task.Titles) == ep;
    sec = Data(:, ...
        Task.Start(ep_ix):(Task.Start(ep_ix)+Task.Frames(ep_ix)-1));
    Segment = [Segment, sec];
end

% If buffer frames are wanted
if isempty(Post_Frames) == false
    ep_ix = string(Task.Titles) == Epoch(end);
    sec = Data(:, ...
        (Task.Start(ep_ix)+Task.Frames(ep_ix)-1):...
        (Task.Start(ep_ix)+Task.Frames(ep_ix)+Post_Frames));
    Segment = [Segment, sec];
end

% Generating table output
Output = table();

% Populating the output table
vars = string(SourceTable.Properties.VariableNames);
vars = vars(vars ~= Variable);
for v = vars
    Output.(v) = SourceTable.(v);
end

% Saving the segment
title = Variable;% strcat(join(Episode, ' + '), " (", Variable, ")");
Output.(title) = Segment;
end

