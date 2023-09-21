%% LOAD THE REQUIRED DATASETS
% 1 - Adding the project path to the searchpath
% Experiment properties
output_path = "C:\Users\lpm97\OneDrive\Documentos\Documentos\Laboratorio de neuro\Git paper\output AUC binn";
% Adding functions to search path
addpath(genpath(output_path))
% 2 - Creating a new folder for the query
target_path = output_path + "\GlobalFluorescence Binning " + ...
    string(datetime(floor(now),'ConvertFrom','datenum'));
mkdir(target_path)
load("ExperimentData.mat")
Experiment.Project.Outputpath=output_path;
%% 1 - Global fluorescence changes
% Extracting global flu
    % Setting the parameters for the specific query
    Data = Experiment.FC;
    Task = Experiment.FC.Task;
    
    % Extracting the palette and sexes
    GroupBy = string(Experiment.Project.Sexes);
    Palette = Experiment.Project.Palette;
    
    % Running the function
    Fluorescence = F_PopulationFluorescence(Data, GroupBy, "Sex", ...
        output_path);
    exportgraphics(gcf, strcat(target_path, "\Outliers.pdf"),  ...
        'ContentType','vector')
% Computing from mean fluorescence DeltaF/F0 - Shifting the traces
    
    % Parameters for the function
    RefEpoch = [];
    Method = "Pctile 50"; % "Mean", "Pctile n", "None"
    Scaling = true;
    
    % Actual shifting the traces
    Fluorescence = F_ShiftTrace(Fluorescence, 'Flu', [], ...
        Method, Task, Scaling);
%% 2 - Binning
% Actually binning
    % Parametters for the binning
    InputTable = Fluorescence;
    Variable = 'Flu';
    BinWidth = 5; %  In seconds
    % Actual binning
    [Binned, Repeated] = F_BinTraces(InputTable, Variable, BinWidth, ...
        Experiment);
        % Where binned are the mean global trace value for the specific bin
        % Frame-by-frame value of their corresponding bin
        % If we have 4-frame-long bins:
            % Binned = [1, 2, 3, 4]
            % Repeated = [1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4]
    % Saving
    Fluorescence.Binned = Binned;
    Fluorescence.RepBinned = Repeated;
%%
% Preallocate space for regression coefficients
n_animals = size(Binned, 1);
slopes = zeros(n_animals, 1);
intercepts = zeros(n_animals, 1);
% For each animal, perform a linear regression on binned data
for i = 1:n_animals
    % Define x and y for regression
    x = (BinWidth/2:BinWidth:(BinWidth/2 + BinWidth*(size(Binned, 2)-1)))'; % x-axis is the midpoint of each bin
    y = Binned(i, :)';
    
    % Perform regression
    p = polyfit(x, y, 1); % p(1) is slope, p(2) is intercept
    
    % Save results
    slopes(i) = p(1);
    intercepts(i) = p(2);
end
% Saving
Fluorescence.Slopes = slopes;
Fluorescence.Intercepts = intercepts;
%%
% Parámetros iniciales
BinWidth = 5; % En segundos (o el valor que hayas definido previamente)
BinWidth_Frames = Experiment.Project.RF * BinWidth;
% Extraer nombres de los animales
animalNames = fieldnames(Experiment.FC);
% Inicializar una estructura para guardar los resultados
PeakBinned = struct();
for i = 1:length(animalNames)-1  % -1 para excluir la variable Task
    animalName = animalNames{i};
    
    % Extraer datos de Filt para el animal actual
    data = Experiment.FC.(animalName).Filt;
    
    % Identificar picos usando islocalmax
    peaks = islocalmax(data, 2); % Asumiendo que quieres identificar máximos a lo largo de las columnas
    
    % Inicializar una matriz para guardar el número total de picos por bin
    totalPeaksPerBin = zeros(1, floor(size(data, 2) / BinWidth_Frames));
    
    % Iterar sobre cada bin y calcular el total de picos
    for j = 1:floor(size(data, 2) / BinWidth_Frames)
        start_ = (j-1)*BinWidth_Frames + 1;
        end_ = j*BinWidth_Frames;
        
        % Extraer datos para el bin actual
        binData = peaks(:, start_:end_);
        
        % Contar el número total de picos por neurona en el bin actual
        totalPeaksPerBin(j) = sum(sum(binData));
    end
    
    % Guardar el promedio de picos por neurona en la estructura PeakBinned
    PeakBinned.(animalName) = totalPeaksPerBin / size(data, 1);  % Divide por el número de neuronas
end
%% 
% Parámetros iniciales para el eje x (basado en la longitud de los bins)
x = (BinWidth/2:BinWidth:(BinWidth/2 + BinWidth*(size(PeakBinned.M2, 2)-1)))'; % x-axis es el punto medio de cada bin
% Extraer nombres de los animales
animalNames = fieldnames(PeakBinned);
% Inicializar una estructura para guardar los coeficientes del ajuste lineal
FitCoefficients = struct();
for i = 1:length(animalNames)
    animalName = animalNames{i};
    
    % Extraer el número de picos por bin para el animal actual
    y = PeakBinned.(animalName);
    
    % Realizar el ajuste lineal
    coefficients = polyfit(x, y, 1); % Grado 1 para ajuste lineal
    
    % Guardar los coeficientes en la estructura FitCoefficients
    FitCoefficients.(animalName) = coefficients;
end