function [Out] = F_RunRScript(script)

threshold = 'FALSE'; %Needed to run correctly the R code 
command = sprintf('C:/PROGRA~1/R/R-43~1.0/bin/Rscript "%s" "%s" %s', script, threshold) 
system(command)

end 
