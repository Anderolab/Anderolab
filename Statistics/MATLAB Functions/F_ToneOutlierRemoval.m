function[Output] = F_ToneOutlierRemovalPlusSTATS(DataTable,Groups,Variable, ID)

    DataTable2 = readtable(DataTable);
    DataTable2.(ID) = string(DataTable2.(ID));
    DataTable2.(Groups) = string(DataTable2.(Groups));
    
    % Tables.one = DataTable2(:,[1 2 3]);
    % Tables.two = DataTable2(:,[1 2 4]);
    % Tables.three = DataTable2(:,[1 2 5]);
    % Tables.four = DataTable2(:,[1 2 6]);
    % Tables.five = DataTable2(:,[1 2 7]);

    for x = 3:width(DataTable2)
        Tables{x-2} = DataTable2(:,[1 2 x])
    end 
    
    % Variables For Grubbs Test
    GroupBy = [Groups];
    Variable1 = [Variable];
    ScriptPath = strcat('', pwd, '\GrubbsTestOutliersMock.R'); % R script path

    threshold = 'FALSE'; % Needed to run correctly the R code
    for i = 1:width(Tables)
        writetable(Tables{1,i}, "TableTunes.csv")
        TableTunes = strcat('', pwd, '\TableTunes.csv');
        command = ['"' 'C:/PROGRA~1/R/R-43~1.0/bin/Rscript" "' ScriptPath ...
        '" "' TableTunes '" "' GroupBy '" "' Variable1 '" "' threshold '"'];
        system(command); 
        my_field = strcat('tone',num2str(i));
        Outliers.(my_field) = readtable("OutliarsTable.csv");
    end
    
    fn = fieldnames(Outliers);
    for i = 1:length(fn)
        if width(Outliers.(fn{i})) ~= 1 
            Output = Outliers.(fn{i});
            for x = height(Output):-1:1;
                outlier = string(table2array(Output(x,2)));
                [index,~] = find(DataTable2{:,:}==outlier);
                DataTable2(index,:) = []; % deletion
            end
        end     
    end
    
    writetable(DataTable2, "DataOutliersRemovedByTone.csv");
    
    % Variables For Statistical Analysis
    GroupBy = [Groups];
    Variable1 = [Variable];
    Id = [ID];
    DataGlobalOutliers = strcat('', pwd, '\DataOutliersRemovedByTone.csv');
    NewFolderPath = strcat('', pwd, '\StatResultsOutliersByTone');
    ScriptPath = strcat('', pwd, '\BoxTestAndANOVA_FC.R'); % R script path
    
    % Run R Script of Analysis
    threshold = 'FALSE'; % Needed to run correctly the R code
    command = ['"' 'C:/PROGRA~1/R/R-43~1.0/bin/Rscript" "' ScriptPath ...
        '" "' DataGlobalOutliers '" "' GroupBy '" "' Variable1 '" "' Id '" "' NewFolderPath ...
        '" "' threshold '"'];
    system(command);
    
    movefile("DataOutliersRemovedByTone.csv", NewFolderPath)
    delete("OutliarsTable.csv", "TableTunes.csv", "Rplots.pdf")

end 