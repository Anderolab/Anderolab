function [Out] = F_MakeCSVForStat(GroupBy, Variable1, Variable2, ...
    Id, DataSetPath, NewFolderPath, DatasetPath_ChiTest)

%Makes matrix and output CSV file 
m = ["Group by:" GroupBy;"First Variable: " Variable1;"Second Variable: " Variable2; 
    "Name od ID column: " Id; 
    "Path for the Dataset: " DataSetPath; 
    "Path for the new folder: " NewFolderPath,
    "Path for ChiTest Dataset: " DatasetPath_ChiTest] 

writematrix(m,'VariablesForStatAnalysis.csv') 

end 