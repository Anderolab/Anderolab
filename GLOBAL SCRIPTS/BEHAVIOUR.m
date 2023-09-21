WorkEnvBH
close all
%%
close all
load("batlow.mat")
subplot(1, 4, [1, 2])
yl = [];
Session = 'FC';
Freez = F_PercentageFreez(Experiment.(Session).Behaviour, false);

% Grouping
Epochs = string(Freez.Properties.VariableNames);
EOI = ['Hab', Epochs(contains(Epochs, 'Tone'))];% Epoch of interest
F_GroupBehaviour(Experiment, Freez, EOI, 'Sexes', Session)
title(Session)
ylabel("% freezing")
%%
subplot(1, 4, 3)
Session = 'FE1';
Freez = F_PercentageFreez(Experiment.(Session).Behaviour, false);

% Grouping
Epochs = string(Freez.Properties.VariableNames);
EOI = Epochs(contains(Epochs, 'Tone'));% Epoch of interest
Groups = {EOI(1:5); EOI(6:10); EOI(10:end)};
GroupName = ["Tones 1-5", "Tones 6-10", "Tones 11-15"];
Freez = F_MergeTimepoints(Freez, Groups, GroupName);
F_GroupBehaviour(Experiment, Freez, ['Hab', GroupName], 'Sexes', Session)
title(Session)


subplot(1, 4, 4)
Session = 'FE2';
Freez = F_PercentageFreez(Experiment.(Session).Behaviour, false);

% Grouping
Epochs = string(Freez.Properties.VariableNames);
EOI = Epochs(contains(Epochs, 'Tone'));% Epoch of interest
Groups = {EOI(1:5); EOI(6:10); EOI(10:end)};
GroupName = ["Tones 1-5", "Tones 6-10", "Tones 11-15"];
Freez = F_MergeTimepoints(Freez, Groups, GroupName);
F_GroupBehaviour(Experiment, Freez, ['Hab', GroupName], 'Sexes', Session)
title(Session)
linkaxes(findall(gcf,'type','axes'), 'y');
f_ = gcf;
f_.Position = [680, 400, 640, 490]
exportgraphics(gcf, 'Behaviour.pdf', "ContentType", 'vector')

%%
close all
colormap(batlow)
subplot(3, 1, 1)
Session = 'FC';
Freez = F_PercentageFreez(Experiment.(Session).Behaviour, true);

subplot(3, 1, 2)
Session = 'FE1';
Freez = F_PercentageFreez(Experiment.(Session).Behaviour, true);

subplot(3, 1, 3)
Session = 'FE2';
Freez = F_PercentageFreez(Experiment.(Session).Behaviour, true);
exportgraphics(gcf, 'DDD.pdf', "ContentType", 'vector')
%%
% Merging epochs



F_GroupBehaviour(Experiment, Freez, ['Hab', GroupName], 'Sexes', Session)

