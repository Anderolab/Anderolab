function [dist_pre_con,dist_tone_con,dist_pen_con,Output, Responses] = F_GetTunningProps(Experiment, Task, TunningEpochs, ...
    ReferenceEpochs, Dataset, Iterations)
%F_GETTUNNINGPROPS Summary of this function goes here
%   Detailed explanation goes here
%% STEP 1 - GENERATING THE STORAGE OUTPUTS
ExpectedTone = []; % From shuffling
ExpectedPretone=[];
ObservedTone = []; % From original
ObservedPretone=[];
StandDevsTone = []; % From shuffling
StandDevsPretone=[];
STD_DistanceTone = []; % Distance in STDs
STD_DistancePretone = [];
STD_DistancePenal=[];
Active = []; % Active neuron traces (original)
Inhibited = []; % Inhibited neuron traces (original)
AllTraces = []; % All traces (original)
Report = []; % Reporting the user of the progress
lengths = []; % Lenths for all sessions
ResponseType=[];
AnimalPerNeuron=[];
SexPerNeuron=[];



OutputPath = Experiment.Project.Outputpath; % Gathering the output path
                        % Saved in the experiment struct
save_path = strcat(OutputPath, "\", join(TunningEpochs(1,:), " - "), ...
    " Tunning ", string(datetime(floor(now),'ConvertFrom','datenum'))); %#ok<TNOW1> 
prompt = strcat("   Results will be saved @ ", save_path);
Report = [Report; prompt; ""];
fprintf('%s\n', prompt)
disp([save_path])
% And the storage site
mkdir(save_path)


%% STEP 2 - IDENTIFYING OUTLIERS BEFORE PERFORMING THE TEST
% Identifying the animals
an_s = fieldnames(Dataset);
Animals = string(regexp(string(an_s(1:end-1)), '\d*', 'Match'));

% Finding all trial lengths;
c = 1;
for animal = Animals.'
    lengths(c) = size(Dataset.('M'+animal).("Raw"), 2);
    c = c+1;
end

% Identifying outliers
outliers = isoutlier(lengths, "mean");

% Reporting outliers to the user
for i = 1:sum(outliers)
    out_an = Animals(outliers);
    out_len = string(lengths(outliers));
    prompt = strcat("       Animal ", ...
        out_an(i), " was identified as an outlier with ",  ...
        out_len(i), ' frames.');
    fprintf('%s\n', prompt)
    Report = [Report; "   OUTLIER DETECTION:"; prompt];
end

% Visualising the outliers
boxplot(lengths)
xticks([])
ylabel("Frames (n)")

% Croppig all sessions and notifying the user
len = min(lengths);
prompt = strcat("       All sessions will be cropped to ", ...
    num2str(len), " frames.");
Report = [Report; prompt; ""];
fprintf('%s\n', prompt)

% Removing outliers
Animals = Animals(outliers == 0);

% Saving current figure
savename = strcat(save_path, "\Tunning - Outliers.pdf");
exportgraphics(gcf, savename, "ContentType", "vector")

%% STEP 3 - GENERATING MAIN USER'S AND STATISTICS OUTPUT TABLE
% Generating the output table
Output = table();
Output.Animal = double(Animals);
Sexes = string(Experiment.Project.Sexes).';
Output.Sex = Sexes(Output.Animal);

% And the variable names for the table storage
if length(TunningEpochs)>1
    Excited_ColName = strcat(join(TunningEpochs(1,:), " & "), " Excited IX");
    Inhibited_ColName = strcat(join(TunningEpochs(1,:), " & "), ...
        " Inhibited IX");
    ExcitedP_ColName = strcat("Probability of ", ...
        join(TunningEpochs(1,:), " & "), " Excited_1");
    InhibitedP_ColName = strcat("Probability of ", ...
        join(TunningEpochs(1,:), " & "), " Inhibited_2");
else 
    Excited_ColName = strcat(TunningEpochs(1,:), " Excited");
    Inhibited_ColName = strcat(TunningEpochs(1,:)," Inhibited");
    ExcitedP_ColName = strcat("Probability of ", TunningEpochs(1,:), ...
        " Excited_1");
    InhibitedP_ColName = strcat("Probability of ", TunningEpochs(1,:), ...
        " Inhibited_2");
end

% Generating the storage variable columns
Output.(Excited_ColName) = repelem({""}, length(Animals)).'; %#ok<STRSCALR> 
Output.(Inhibited_ColName) = repelem({""}, length(Animals)).'; %#ok<STRSCALR> 
Output.(ExcitedP_ColName) = zeros(length(Animals), 1);
Output.(InhibitedP_ColName) = zeros(length(Animals), 1);

% For figures
Ep_ = {ExcitedP_ColName, InhibitedP_ColName};
%% STEP 4 - GENERATING INDEXES FOR FRAMES OF INTEREST
% Identifying frames corresponding to the time of interest
TOITone = []; % Time of interest
TOIPre = [];
for Epoch = TunningEpochs(1, :) % What? - % TunningEpochs(2,:)
    ixTone = string(Task.Titles) == Epoch;
    TOITone = [TOITone, Task.Start(ixTone):(Task.Start(ixTone)+Task.Frames(ixTone)-30)];
end
for Epoch = TunningEpochs(1, :) % What? - % TunningEpochs(2,:)
    ixPre = find(string(Task.Titles) == Epoch, 1)-1;
    TOIPre = [TOIPre, (Task.Start(ixPre)+Task.Frames(ixPre)-30*15):(Task.Start(ixPre)+Task.Frames(ixPre))];
end

% Identifying frames corresponding to the frames of reference
TOR = []; % Frames of reference
for Epoch = ReferenceEpochs
    ix = string(Task.Titles) == Epoch;
    TOR = [TOR, Task.Start(ix):(Task.Start(ix)+Task.Frames(ix))];
end

% Generating the binary for the comparison
BinarisedTask = repelem("0", len);
BinarisedTask(TOITone) = "1";
BinarisedTask(TOIPre)="2";
area(str2double(BinarisedTask)*2, "FaceAlpha", .3)
prompt = "   Threshold will be of 1.96 SDs from mean (95% CI)";
Report = [Report; ""; prompt];
fprintf('%s\n', prompt)
%% STEP 5 - PERFORMING THE TEST
% Looping through animals
c = 1; % Counter function
for an = Animals.'
    wb_ = waitbar(0, ...
        strcat("Identifying responsive neurons in animal ", num2str(an)));
    prompt = strcat("   Processing animal ", num2str(an));
    fprintf('%s\n', prompt)
    Report = [Report; prompt];
    % Gathering the animal specific data
    Maxims = double(islocalmax(Dataset.(strcat('M', ...
        num2str(an))).Filt(:, 1:len), 2));

    % subplot(2,2,1)
    % if str2double (an)==2
    %     a=sum(Maxims == 1, 2);
    % 
    %     [~, indice] = max(a);
    %     disp(indice)
    %     area(str2double(BinarisedTask)*2, "FaceAlpha", .3)
    %     hold on
    %     imagesc(Maxims(indice,:))
    % end
    Peaks = string(Maxims);
    Intersect = BinarisedTask+Peaks;
    % Identifying the parameters for each neuron 
    oo = sum(Intersect == "00", 2);
    lo = sum(Intersect == "10", 2);
    ol = sum(Intersect == "01", 2);
    ll = sum(Intersect == "11", 2);
    sl=sum(Intersect == "21", 2);
    so=sum(Intersect == "20", 2);
    mod_oo=oo+so;
    mod_ll=ll;
    mod_ol=ol+sl;
    mod_lo=lo;
    
    
    ObvTone = F_ComputePhi(oo+so, ol+sl, lo, ll);
    ObvPretone = F_ComputePhi(oo, ol, so, sl);
    % Computing Phi
    ObservedTone = [ObservedTone; ObvTone];
    ObservedPretone=[ObservedPretone; ObvPretone];

    % Rescaling the traces for the randomisation
    MaximsTone = Maxims;
    MaximsPretone = Maxims;

    % Eliminar los periodos de interés opuestos
    %MaximsTone(:, TOIPre) = [];
    MaximsPretone(:, TOITone) = [];

    lenTone = size(MaximsTone, 2);
    lenPretone = size(MaximsPretone, 2);
    

    % Iterating to attain the expected
    Expected_ScoresTone = gpuArray(zeros(size(MaximsTone, 1), Iterations));
    Expected_ScoresPretone = gpuArray(zeros(size(MaximsPretone, 1), Iterations));

    PeaksTone=string(MaximsTone);
    PeaksPretone=string(MaximsPretone);
    % Generate BinarisedTask for Tone
    BinarisedTaskTone = BinarisedTask;
    %BinarisedTaskTone(TOIPre) = [];
    
    % Generate BinarisedTask for Pretone
    BinarisedTaskPretone = BinarisedTask;
    BinarisedTaskPretone(TOITone) = [];


    for iter = 1:Iterations

        waitbar(iter/Iterations);
        Rand_PeaksTone = PeaksTone(:, randperm(lenTone));
        
        Rand_PeaksPretone = PeaksPretone(:, randperm(lenPretone));
        % if iter==1 & str2double(an)==2
        %     subplot(2,2,2)
        %     area(str2double(BinarisedTaskTone)*2, "FaceAlpha", .3)
        %     hold on
        %     imagesc(str2double(Rand_PeaksTone(indice,:)))
        %     subplot(2,2,3)
        %     area(str2double(BinarisedTaskPretone)*2, "FaceAlpha", .3)
        %     hold on
        %     imagesc(str2double (Rand_PeaksPretone(indice,:)))
        % end
        Intersect_RandTone = BinarisedTaskTone + Rand_PeaksTone;
        Intersect_RandPretone = BinarisedTaskPretone + Rand_PeaksPretone;
    
        % For Tone
        oo_tone = sum(Intersect_RandTone == "00", 2);
        lo_tone = sum(Intersect_RandTone == "10", 2);
        ol_tone = sum(Intersect_RandTone == "01", 2);
        ll_tone = sum(Intersect_RandTone == "11", 2);
        sl_tone=sum(Intersect_RandTone == "21", 2);
        so_tone=sum(Intersect_RandTone == "20", 2);
        Expected_ScoresTone(:, iter) = F_ComputePhi(oo_tone+so_tone, ol_tone+sl_tone, lo_tone, ll_tone);
    
        % For Pretone
        oo_pre = sum(Intersect_RandPretone == "00", 2);
        so_pre = sum(Intersect_RandPretone == "20", 2);
        ol_pre = sum(Intersect_RandPretone == "01", 2);
        sl_pre = sum(Intersect_RandPretone == "21", 2);
        Expected_ScoresPretone(:, iter) = F_ComputePhi(oo_pre, ol_pre, so_pre, sl_pre);

    end
    ExpTone = mean(Expected_ScoresTone, 2);
    ExpPre = mean(Expected_ScoresPretone, 2);
    ExpectedTone = [ExpectedTone; ExpTone];
    ExpectedPretone = [ExpectedPretone; ExpPre];
    SDsTone = std(Expected_ScoresTone, [], 2);
    SDsPretone = std(Expected_ScoresPretone, [], 2);
    StandDevsTone = [StandDevsTone; SDsTone];
    StandDevsPretone = [StandDevsPretone; SDsPretone];
    DistancesTone = (ObvTone-ExpTone)./SDsTone;
    disp(size(DistancesTone));
    DistancesPretone = (ObvPretone-ExpPre)./SDsPretone;
    STD_DistanceTone = [STD_DistanceTone; DistancesTone];
    STD_DistancePretone = [STD_DistancePretone; DistancesPretone];
    %yuju=DistancesPretone./DistancesTone;
    Distances_Penal = DistancesTone;
    condition = (sign(DistancesPretone) == sign(DistancesTone) & (DistancesTone > DistancesPretone)& abs(DistancesPretone)>1.96 & abs(DistancesTone)>1.96);

    Distances_Penal(condition)=(1-DistancesPretone(condition)./DistancesTone(condition)).*DistancesTone(condition);
    if isinf(Distances_Penal(condition))
        "PRETONE"
        DistancesPretone(condition)
        
        "TONE"
        DistancesTone(condition)
    end
    dist_pen_con=Distances_Penal(condition)
    dist_pre_con=DistancesPretone(condition)
    dist_tone_con=DistancesTone(condition)
    % for d=1:length(DistancesTone)
    %     if DistancesPretone(d)>0 && DistancesTone(d)>0 && DistancesTone(d)>DistancesPretone(d)
    %         Distances_Penal_tra(d)=(1-DistancesPretone(d)./DistancesTone(d)).*DistancesTone(d);
    %     else
    %         Distances_Penal_tra(d)=DistancesTone(d);
    % 
    %     end
    % 
    % 
    % end
    % Distances_Penal=Distances_Penal_tra';
    % culo=size(Distances_Penal)
    STD_DistancePenal=[STD_DistancePenal;Distances_Penal];

    % Saving the activated and inhibited neurons
    Active = [Active; ...
        Dataset.(strcat('M', num2str(an))).Filt(Distances_Penal > 1.96, 1:len)];
    Inhibited = [Inhibited; ...
        Dataset.(strcat('M', num2str(an))).Filt(Distances_Penal < -1.96, 1:len)];
    % All traces
    AllTraces = [AllTraces; ...
        Dataset.(strcat('M', num2str(an))).Filt(:, 1:len)];

    % Saving
    Excit = find(Distances_Penal > 1.96);
    Inhibit = find(Distances_Penal < -1.96);
    Output.(Excited_ColName)(c) = {Excit};
    Output.(Inhibited_ColName)(c) = {Inhibit};
    Output.(ExcitedP_ColName)(c) = 100*length(Excit)/size(Maxims, 1);
    Output.(InhibitedP_ColName)(c) = 100*length(Inhibit)/size(Maxims, 1);
    
        % For the frequency test
        Tunning_ = repelem("Unresponsive", length(Distances_Penal));
        Tunning_(Excit) = "Excited";
        Tunning_(Inhibit) = "Inhibited";
        ResponseType = [ResponseType, Tunning_];

        % Saving the animal
        AnimalPerNeuron = [AnimalPerNeuron, ...
            repelem(an, length(Distances_Penal))];
        SexPerNeuron = [SexPerNeuron, ...
            repelem(Experiment.Project.Sexes{double(an)}, ...
            length(Distances_Penal))];

    prompt = strcat("       ",  num2str(length(Excit)), ...
        " stimulus-excited neurons have been identified for animal ", ...
        num2str(an));
    fprintf('%s\n', prompt)
    Report = [Report; prompt];
    prompt = strcat("       ",  num2str(length(Inhibit)), ...
        " stimulus-inhibited neurons have been identified for animal ", ...
        num2str(an));
    fprintf('%s\n', prompt)
    Report = [Report; prompt];    
    c = c +1;
    close(wb_);

end

%% STEP 6 - GENERATING THE VISUALISATIONS
% First figure - Methods
    [~, sort_ix] = sort(ObservedTone);
    F_FillArea(ExpectedTone(sort_ix).', (StandDevsTone(sort_ix).*1.96).', ...
        'k', 1:length(ExpectedTone(sort_ix)))
    hold on
    plot(ExpectedTone(sort_ix), "Color", 'k')
    hold on
    plot(ObservedTone(sort_ix), "Color", 'r', "LineWidth", 2)
    O = ObservedTone(sort_ix);
    STD_Sorted = STD_DistanceTone(sort_ix);
    sig_ix = find(abs(STD_Sorted) > 1.96);
    scatter(sig_ix, O(sig_ix), 20, 'K', 'filled')
    hold off
    legend(["95% CI", "Expected", "Observed", "Significant"], ...
        "Location","northwest");
    xlim([1, length(ObservedTone)])
    xlabel("Neuron", "FontSize", 12)
    ylabel("\phi Coefficient", "FontSize", 12)
    set(gcf,'Position',[400 100 300 400])
    box off
    hold off
    savename = strcat(save_path, "\Tunning - Neuron identification.pdf");
    exportgraphics(gcf, savename, "ContentType","vector")

    
% Second figure - Sample
    % For the active neurons - Calculating the error
    close all
    sem_ = std(Active, [], 1)./sqrt(size(Active, 1));
    % Getting the axes
    F_FillArea(mean(Active, 1), sem_, [212, 100, 66]./255, ...
        1:length(Active))
    yl = ylim();
    
    % Labelling areas of interest
    area(double(BinarisedTask), 'EdgeColor','none', "FaceAlpha", .3, ...
        'FaceColor', 'k')
    hold on

    % Visualising the errors for the active neurons
    F_FillArea(mean(Active, 1), sem_, [212, 100, 66]./255, ...
        1:length(Active))
    % Viewing active neurons
    plot(mean(Active, 1), "Color", [212, 100, 66]./255)
    hold on
    F_FillArea(mean(Inhibited, 1), sem_, [76, 148, 199]./255, ...
        1:length(Inhibited))
    plot(mean(Inhibited, 1), "Color", [76, 148, 199]./255)
    
    ylim([0, max(yl)])
    xlim([1, length(Active)])
    set(gcf, 'Position', [400, 100, 1200, 500])
    box off
    xlabel("Time (Frames)")
    ylabel("Mean Fluorescence")
    legend([join(TunningEpochs(1,:), ' & '), "", "Excited", "", "Inhibited"])
    hold off
    savename = strcat(save_path, "\Tunning - Group mean.pdf");
    exportgraphics(gcf, savename, "Resolution", 1200)

% Third figure - Individual-neuron level visualisation
    close all
    % Selecting the top five neurons;
    Sorted = sort(STD_DistancePenal);
    TopInactive = AllTraces(STD_DistancePenal < Sorted(21), :);

    % Getting the axes
    area(double(BinarisedTask).*13, 'EdgeColor','none', ...
        "FaceAlpha", .3, 'FaceColor', 'k')
    hold on
    for i = 1:5
        plot(TopInactive(i, :) + i, "Color", [76, 148, 199]./255)
        hold on
    end
    Sorted = flip(Sorted);
    TopActive = AllTraces(STD_DistancePenal > Sorted(21), :);
    for i = 7:11

        plot(TopActive(i-6, :)+i, "Color", [212, 100, 66]./255)
        hold on
    end
    hold off
    ylim([0, 13])
    yticks([3, 9])
    yticklabels({"Inhibited", "Excited"}); %#ok<CLARRSTR> 
    ytickangle(90)
    xlabel("Time (Frames)")
    ylabel("Filtered fluorescence")
    xlim([1, size(AllTraces, 2)])
    box off
    Fig_ = gcf;
    Fig_.Position = [400, 100, 600, 500];
    savename = strcat(save_path, "\Tunning - Sample Neurons.pdf");
    exportgraphics(gcf, savename, "ContentType", "vector")
    
    sample_neurons = [TopActive(1:10, TOIPre(1):[TOITone(2)+60]); ...
        TopInactive(1:10, TOIPre(1):[TOITone(2)+60])];
    imagesc(sample_neurons)

% Fourth figure - Ratios
    % Identifying the percentage of neuron tunning
    "FFFFf"
    close all
    Properties = [sum(Distances_Penal > 1.96, 'All'), ...
        sum(Distances_Penal < -1.96, 'All')];
    Properties(end+1) = length(Distances_Penal) - sum(Properties);
    pie(Properties)
    ax = gca();
    ax.Colormap = [212, 100, 66; 76, 148, 199; 200, 200, 200]./255; 
    legend(["Excited", "Inhibited", "Irresponsive"], "Location", ...
        "northeast")
    F_ = gcf;
    F_.Position = [400, 100, 390, 320];
    savename = strcat(save_path, "\Tunning - Ratios.pdf");
    exportgraphics(gcf, savename, "ContentType", "vector")
 
% Fifth figure - Distinction between conditions
    % Generating the first subplot
    close all
    yls = [];
    subplot(1, 2, 1)

    boxplot(Output.(ExcitedP_ColName), Output.Sex, ...
        "Colors", 'k');
    colours = ...
        flip(cell2mat(Experiment.Project.Palette(unique(Output.Sex, ...
        'stable'))), 1);
    f = findobj(gca, 'Tag', 'Box');
    for i = 1:length(f)
        patch(get(f(i),'XData'),get(f(i),'YData'),colours(i,:), ...
            'FaceAlpha',.6);
    end
    h=findobj('LineStyle','--'); set(h, 'LineStyle','-');
    box off
    hold on
    ylabel("Tunned neurons (%)")
    xlabel("Excited")
    yls = [yls; ylim()];

    subplot(1, 2, 2)
    boxplot(Output.(InhibitedP_ColName), Output.Sex, ...
        "Colors", 'k');
    f = findobj(gca, 'Tag', 'Box');
    for i = 1:length(f)
        patch(get(f(i),'XData'),get(f(i),'YData'),colours(i,:), ...
            'FaceAlpha',.6);
    end
    h=findobj('LineStyle','--'); set(h, 'LineStyle','-');
    xlabel("Inhibited")
    h = gca;
    h.YAxis.Visible = 'off';
    box off
    yls = [yls; ylim()];
    ylim([min(yls(:, 1)), max(yls(:, 2))])

    
    savename = strcat(save_path, "\Tunning - Group Differences.pdf");
    exportgraphics(gcf, savename, "ContentType", "vector")
    close all

%% STEP 7 - CLOSING
% Saving the report
savename = strcat(save_path, "\Report.txt");
writelines(Report,savename);
% And the data
savename = strcat(save_path, "\TunningOutput.mat");
Tunning = Output;
save(savename, "Tunning");

% Generating the frequency table
Responses = table();
size(SexPerNeuron)
size(ResponseType)
size(AnimalPerNeuron)
Responses.Sex = SexPerNeuron.';
Responses.Animal = AnimalPerNeuron.';
Responses.ResponseType = ResponseType.';

% And saving it
savename = strcat(save_path, "\SingleNeuronResponses.mat");
save(savename, "Responses");
savename_chi = strcat(save_path, "\SingleNeuronResponses.csv");
writetable(Responses, savename_chi);

%% STEP 8 - STATS
% Save Table of Data
% savename = strcat(save_path, "\TunningOutput.csv");
% StatsTable = table();
% StatsTable.Sex = Output.Sex;
% StatsTable.Animal = Output.Animal;
% StatsTable.(ExcitedP_ColName) = Output.(ExcitedP_ColName);
% StatsTable.(InhibitedP_ColName) = Output.(InhibitedP_ColName);
% writetable(StatsTable, savename);
% 
% % Make Table of Variables 
% GroupBy = "Sex";
% Variable1 = ExcitedP_ColName;
% Variable2 = InhibitedP_ColName;
% Id = "Animal";
% DatasetPath = string(savename);
% NewFolderName = "StatResults";
% NewFolderPath = strcat(save_path, "\", NewFolderName);
% DatasetPath_ChiTest = string(savename_chi)
% 
% F_MakeCSVForStat(GroupBy, Variable1, Variable2, Id, DatasetPath, NewFolderPath, DatasetPath_ChiTest);
% 
% % Run R Script of Analysis
% ScriptPath = strcat('', pwd, '\BoxTestAndANOVA.R')
% F_RunRScript(ScriptPath)


end

