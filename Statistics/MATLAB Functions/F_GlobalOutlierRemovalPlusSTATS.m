function[Output] = F_GlobalOutlierRemovalPlusSTATS(DataTable,Groups,Variable, ID)

    % Variables For Grubbs Test
    GroupBy = [Groups];
    Variable1 = [Variable];
    ScriptPath = strcat('', pwd, '\GrubbsTestOutliersMock.R'); % R script path
    
    % Run R Script of Analysis
    threshold = 'FALSE'; % Needed to run correctly the R code
    command = ['"' 'C:/PROGRA~1/R/R-43~1.0/bin/Rscript" "' ScriptPath ...
        '" "' DataTable '" "' GroupBy '" "' Variable1 '" "' threshold '"'];
    system(command);
    
    % Outliers Table 
    Output = readtable("OutliarsTable.csv");
    
    DataTable2 = readtable(DataTable);
    DataTable2.(ID) = string(DataTable2.(ID));
    DataTable2.(Groups) = string(DataTable2.(Groups));
    
    % To Remove Outliers 
    if width(Output) ~= 1;
        for x = height(Output):-1:1;
            outlier = string(table2array(Output(x,2)));
            [index,~] = find(DataTable2{:,:}==outlier);
            DataTable2(index,:) = []; % deletion
        end 
    end
    
    writetable(DataTable2, "DataGlobalOutliersRemoved.csv");
    
    % Variables For Statistical Analysis
    GroupBy = [Groups];
    Variable1 = [Variable];
    Id = [ID];
    DataGlobalOutliers = strcat('', pwd, '\DataGlobalOutliersRemoved.csv');
    NewFolderPath = strcat('', pwd, '\StatResultsGlobalOutliers');
    ScriptPath = strcat('', pwd, '\BoxTestAndANOVA_FC.R'); % R script path
    
    % Run R Script of Analysis
    threshold = 'FALSE'; % Needed to run correctly the R code
    command = ['"' 'C:/PROGRA~1/R/R-43~1.0/bin/Rscript" "' ScriptPath ...
        '" "' DataGlobalOutliers '" "' GroupBy '" "' Variable1 '" "' Id '" "' NewFolderPath ...
        '" "' threshold '"'];
    system(command);
    
    movefile("DataGlobalOutliersRemoved.csv", NewFolderPath);


end 