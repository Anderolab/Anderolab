MockData = table();
MockData.Sex =["Female"; "Female"; "Female"; "Female"; "Female"; "Female"; "Female"; "Male"; ...
    "Male"; "Male"; "Male"; "Male"; "Male"; "Male"; "Male"; "Male"];
MockData.Values = vertcat(num2cell(randn(4,1)), 30.666, num2cell(randn(5,1)), -40.00, 30.09, num2cell(randn(4,1)));

%%

% Write csv fpr of the data
writecell(MockData.Values, "AbsFlu.csv");

% Run R Script of Analysis
ScriptPath = strcat('', pwd, '\GrubbsTestOutliers.R');
F_RunRScript(ScriptPath);

Output_Global = readtable("OutliarsTable.csv")
subplot(1, 4, [1,2])
boxplot([MockData.Values{:}])
y_ = ylim();

%%

MockData_Fem = MockData(1:7,:);
MockData_Mal = MockData(8:end,:);

%%
% Write csv fpr of the data
writecell(MockData_Fem.Values, "AbsFlu.csv");

% Run R Script of Analysis
ScriptPath = strcat('', pwd, '\GrubbsTestOutliers.R');
F_RunRScript(ScriptPath);

Output_Fem = readtable("OutliarsTable.csv")

subplot(1, 4, 3)
boxplot([MockData_Fem.Values{:}])
ylim(y_)

%%
% Write csv fpr of the data
writecell(MockData_Mal.Values, "AbsFlu.csv");

% Run R Script of Analysis
ScriptPath = strcat('', pwd, '\GrubbsTestOutliers.R');
F_RunRScript(ScriptPath);

Output_Mal = readtable("OutliarsTable.csv")

subplot(1, 4, 4)
boxplot([MockData_Mal.Values{:}])
ylim(y_)




