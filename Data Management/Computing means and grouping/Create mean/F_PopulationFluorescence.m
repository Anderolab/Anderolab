
function[Output] = F_PopulationFluorescence(Data, GroupBy, ...
    GroupingVariableNames, TargetPath)

%% STEP 1 - Alligning all the sessions
% Extracting all available animal numbers
an_s = fieldnames(Data);
Animals = string(regexp(string(an_s(1:end-1)), '\d*', 'Match'));
Report = [];

% Finding all trial lengths
lengths = [];
c = 1;
for animal = Animals.'
    lengths(c) = size(Data.('M'+animal).("Raw"), 2);
    c = c+1;
end
% Identifying outliers via R
[OutlierTable, outliers] = F_GrubbsTest(lengths.');

% Notifying the user and generating the analysis report
for i = 1:sum(outliers)
    out_an = Animals(logical(outliers));
    out_len = string(lengths(logical(outliers)));
    out_p = OutlierTable.P_Values;
    prompt = strcat("    Animal ", ...
        out_an(i), " was identified as an outlier with ",  ...
        out_len(i), " frames. It will be excluded from the rest ", ...
        "of the sample. Grubb's p = ", num2str(out_p(i)));
    fprintf('%s\n', prompt)
    Report = [Report; prompt];
end

% Generating visualisation
boxplot(lengths)
ylabel("Session length (frames)")
xticklabels("All animals")
exportgraphics(gcf, strcat(TargetPath, "\SessionLengthOutliers.pdf"),  ...
    'ContentType','vector')

% Croppig all sessions
len = min(lengths(outliers == 0));

prompt = strcat("    All sessions will be cropped to ", ...
    num2str(len), " frames.");
Report = [Report; prompt; ""];
fprintf('%s\n', prompt)

Animals = Animals(outliers == 0);

%% STEP 2 - Calculating the means
Means = zeros(length(Animals), len);    % Storage
c = 1;  % Counter
for a = Animals.'
    Means(c, :) = mean(Data.('M'+a).('Raw')(:, 1:len));
    c = c+1;
end

%% STEP 3 - Identifying sessions with abnormal fluorescence.
% Computing the mean absolute fluorescence
AbsFlu = mean(abs(Means), 2);

% Perfroming Grubb's via R
[Outs, FluOuts] = F_GrubbsTest(AbsFlu);

% Updating the user and populating the analysis report
if sum(FluOuts) ~= 0
    for line_ = 1:size(Outs, 1)
        ix_ = find(FluOuts, line_);
        an_ = Animals(ix_(end));
        p_ = Outs.P_Values(line_);
        prompt = strcat("    Animal ", an_, " will be excluded from ", ...
            "further analyses since it was identified as an outlier with ", ...
            "an abnormal fluorescence. Grubb's p = ", num2str(p_), ".");
        Report = [Report; prompt];
        fprintf('%s\n', prompt)
    end
end

% Generating and saving a graphic of outliers
boxplot(AbsFlu)
xlim([.5, 1.5])
xticklabels("All animals")
ylabel("Mean absolute fluorescence")
exportgraphics(gcf, strcat(TargetPath, ...
    "\GlobalFluorescenceOutliers.pdf"), 'ContentType','vector')

%% STEP 4 - Generating output
Output = table();
Output.(GroupingVariableNames) = GroupBy(double(Animals)).';
Output.Flu = Means;
Output.Animal = double(Animals);
Output = Output(~FluOuts, :); % Excluding outlier rows.
end