%% Loading the data
% Loading the FC to simulate the first imo (IMO20)
IMOE = load("All_IMO_Data_E.mat");

% Loading the FE1 to simulate the second imo (IMO10)
IMOL = load("All_IMO_Data_L.mat");



%%
% Sex of each animal: - USAR ESTE NO EL Cat.Dat_SEX
Sexes = {"Female", "Female", "Female", "Female", "Female", "Female", ...
    "Female", "Male", "Male", "Male", "Male", "Male", "Male", "Male"}; %#ok<CLARRSTR>

% NOTE - For consistency purposes, let's use this dictionary to settle on
% colours for males and females. Later on (low priority) we can design a
% gui-function that prompts the user to select their own colour scheme.
ColorDict = dictionary(["Male", "Female"], ...
    {[0, 75, 75]./255, [245, 173, 82]./255});

AnimalReg = {IMOE.IMO_Data.AIX, IMOL.IMO_Data.AIX}; % Animal index per neuron
RF = 30; % FPS
EpochNames = {"Early", "Late"}; %#ok<CLARRSTR> 

save_path = "C:\Users\Ander\OneDrive\Documents\MISCELL\TestEnv\Joaquín\testIMOS";
save_path = strcat(save_path, "\Analysis ", string(datetime(floor(now),'ConvertFrom','datenum')));
format long 
 %target_path = output_path + "\GlobalFluorescence FE1" + ...
    %string(datetime(floor(now),'ConvertFrom','datenum'));

mkdir(save_path)

path_stats_folder = strcat(save_path, "\StatResults_");
%% AUC COMPARE

close all
close(gcf)

AUC_Comparison = F_Compare_AUC(IMOE.IMO_Data.FiltTraces, ...
    IMOL.IMO_Data.FiltTraces, RF, Sexes, AnimalReg, EpochNames, true,  ...
    "Mean firing frequency", ColorDict, char(save_path))

image_AUC = strcat(path_stats_folder, "AUC\Sample - AUC Compare.pdf");
exportgraphics(gcf,image_AUC, "ContentType","vector");
report_AUC = strcat(path_stats_folder, "AUC\Sample_AUC_Compare.mat");
save(report_AUC, 'AUC_Comparison');


close(gcf)

%% PEAKS COMPARE
% Testing the function - Output displayed in command.
close all  

Peaks_E=islocalmax(IMOE.IMO_Data.FiltTraces,2);
Peaks_L=islocalmax(IMOL.IMO_Data.FiltTraces,2);

comparisons = {'Epochs_peak_2:Male-Epochs_peak_2:Female', 'Epochs_peak_2:Male-Epochs_peak_1:Male', ...
    'Epochs_peak_2:Female-Epochs_peak_1:Female','Epochs_peak_1:Male-Epochs_peak_1:Female'};

PEAK_Comparison = F_Compare_PEAKS(Peaks_E, ...
    Peaks_L, RF, Sexes, AnimalReg, EpochNames, true, ...
    "Mean firing frequency", ColorDict, char(save_path), comparisons, 'PosthocEpochsGroup""',true)

image_PEAKS = strcat(path_stats_folder, "PEAKS\Sample - Peaks Compare.pdf");
exportgraphics(gcf,image_PEAKS, "ContentType","vector");
report_PEAKS = strcat(path_stats_folder, "PEAKS\Sample_PEAKS_Compare.mat");
save(report_PEAKS, 'PEAK_Comparison');
% close(gcf)

%% Raster visualisation
numframs = 10000;
numNeurons = 50;
  % 1 - Finding the best animal for each sex
    % Storage
    select_animals = zeros(1, length(unique(PEAK_Comparison.Groups)));
    % Iteration
    c = 1;
    for group = unique(PEAK_Comparison.Groups).'
        % Setting up the parameters
        gr_ix = PEAK_Comparison.Groups == group;
        gr_ans = PEAK_Comparison.Animals(gr_ix);
        % Computing DeltaPeaks
        diff = abs(PEAK_Comparison.Epochs_peak(gr_ix, 1) - ...
            PEAK_Comparison.Epochs_peak(gr_ix, 2));
        % Finding best fit
        select_animals(c) = gr_ans(diff == max(diff));
        % Counter funcion
        c = c+1;
    end
    clear c

  % 2 - Generating the rasterplots
    ix = {[1, 3, 5], [2, 4, 6], [11, 13, 15], [12, 14, 16]};
    c = 1;
    hist_max = zeros(1, 4);
    for animal = select_animals
        for contrast = {IMOE, IMOL}
            % Loading the data
            contrast = contrast{1};
            an_dat = ...
                contrast.IMO_Data.FiltTraces(contrast.IMO_Data.AIX == ...
                animal, 1:numframs);
            % Extracting the peaks
            [n, f] = find(islocalmax(an_dat(1:numNeurons, :), 2) == 1);
            subplot(9, 2, ix{c})
            scatter(f./RF, n, 5, 'k', "filled")

            xlim([1, numframs/RF])
            ylim([1, numNeurons])
            xticks([])
            xticklabels([])
            h = gca;
            h.XAxis.Visible = 'off';
            box off
            subplot(9, 2, max(ix{c})+2)
            histogram(f./RF, floor(numframs/(RF*5)), "FaceColor", 'k', ...
                'FaceAlpha', 1)
            xlim([[1, numframs/RF]])
            hist_max(c) = max(ylim());
            c = c+1;
            box off
            
        end
    end
    subplot(9, 2, ix{1})
    title("Early IMO")
    ylabel({char(9792), "Neurons"}, 'FontSize', 11, 'FontWeight', 'bold')
    yticklabels
    subplot(9, 2, ix{2})
    title("Late IMO")
    subplot(9, 2, ix{3})
    ylabel({char(9794), "Neurons"}, 'FontSize', 11, 'FontWeight', 'bold')
    subplot(9, 2, max(ix{3})+2)
    xlabel("Time (s)")
    subplot(9, 2, max(ix{4})+2)
    xlabel("Time (s)")

    % Scaling to equal
    for spt = 1:(c-1)
        subplot(9, 2, max(ix{spt})+2)
        ylim([0, max(hist_max)])
        yticks([0, max(hist_max)])
    end
    
    f_ = gcf;
    f_.Position = [294, 369, 655, 420];

image_Rastr = strcat(path_stats_folder, "PEAKS\Rastr.pdf");
exportgraphics(gcf,image_Rastr,'Resolution',2000);


%% AMPL COMPARE
% Testing the function - Output displayed in command.
close all

AMPL_Comparison = F_Compare_AMPL2(IMOE.IMO_Data.FiltTraces, ...
    IMOL.IMO_Data.FiltTraces, Sexes, AnimalReg, EpochNames, true, ...
    "Mean amplitude frequency", ColorDict,'Raw', char(save_path))

image_AMPL = strcat(path_stats_folder, "AMPL\Sample - Amplitude compare.pdf")
exportgraphics(gcf,image_AMPL, "ContentType", "vector")
report_AMPL = strcat(path_stats_folder, "AMPL\Sample_AMPL_Compare.mat")
save(report_AMPL, 'AMPL_Comparison')

close(gcf)

%%
Compare_AUC_PEAKS(AUC_Comparison, PEAK_Comparison)
%%


%%AUC_PEAKS_AMPL_Comparison=Compare_AUC_PEAKS_AMPL(AUC_Comparison, PEAK_Comparison,AMPL_Comparison)

%% MOVEMENT HIGH ACTIVITY ANALYSIS
close all
% Pre-processing the dataset
% Performing 95%CI rescaling
MovE = F_95CI_Norm(IMOE.IMO_Data.Movement);
%%
MovL = F_95CI_Norm(IMOL.IMO_Data.Movement);

% Concatenating the two structs to generate the model
MaxLen = max(size(MovE, 2), size(MovL, 2));

% Resizing and concatenating
MovE(:, size(MovE, 2):MaxLen) = NaN;
MovL(:, size(MovL, 2):MaxLen) = NaN;
Global = [MovL; MovE];

Model = F_FitGauss(Global, 'gauss2');

% Generating the model and defining the threshold
Model = F_FitGauss(Global, 'gauss2');

% Defining the gaussian fit functions
gauss2 = @(x) Model.a2*exp(-((x-Model.b2)^2)/(2*Model.c2^2));
gauss1 = @(x) Model.a1*exp(-((x-Model.b1)^2)/(2*Model.c1^2));

% Defining the intersect
diff_func = @(x) gauss1(x) - gauss2(x);
LowBound = fzero(diff_func, [Model.b1, Model.b2]);
HighBound = (1.96*Model.c2)/sqrt(2) + Model.b2;
hold on

% Visualising the area of interest
x = -2:.0001:4;
y = zeros(size(x));
y(x > LowBound & x < HighBound) = 1;
area(x, y, 'FaceAlpha', .3);
legnd = string({gca().Legend.String{:}});
legend(["", "", legnd(1:end-1), "High Activity"])

% Setting the X label
mkdir(strcat(path_stats_folder, "MOV"))
image_MOD = strcat(path_stats_folder, "MOV\Movement model.pdf")
xlabel("(Act-M_d)/95CI")
exportgraphics(gcf,image_MOD,"ContentType", "vector")

% Setting the animal register
AnimalReg = {IMOE.IMO_Data.Animals; IMOL.IMO_Data.Animals};


% Comparing frames in high activity
HM_Comparison = F_Compare_MOV(MovE, MovL, Sexes, AnimalReg, EpochNames, ...
    true, "Probability of high activity", ColorDict, LowBound, HighBound,char(save_path))

image_MOV = strcat(path_stats_folder, "MOV\Sample - MovementModel.pdf")
exportgraphics(gcf,image_MOV,"ContentType", "vector")
report_MOV = strcat(path_stats_folder, "MOV\Sample_HM_Compare.mat")
save(report_MOV, "HM_Comparison");


%% Generating the model to identify the high movement epochs
close all
Model = F_FitGauss(Global, 'gauss2');
LowBound = Model.b1+0.67*(Model.c1/sqrt(2));
HighBound = Model.b2-0.67*(Model.c2/sqrt(2));

hold on
xline(Model.b1+0.674*(Model.c1/sqrt(2)), 'Color', 'r', 'LineWidth', 2, ...
    'LineStyle', ':')
xline(Model.b2-0.674*(Model.c2/sqrt(2)), 'Color', 'r', 'LineWidth', 2, ...
    'LineStyle', ':')
hold off

image_EV = strcat(path_stats_folder, "MOV\Sample - EventDefiningModel.pdf")
legnd = string({gca().Legend.String{:}});
legend(["", "", legnd(1:end-1), "LowBound", "HighBound"])

exportgraphics(gcf,image_EV,"ContentType", "vector")
% 
% exportgraphics(gcf,'Sample - EventDefiningModel.pdf',"ContentType", "vector")


%%
Data1 = MovE;
Data2 = MovL;

% Sex of each animal: - USAR ESTE NO EL Cat.Dat_SEX
Sexes = {"Female", "Female", "Female", "Female", "Female", "Female", ...
    "Female", "Male", "Male", "Male", "Male", "Male", "Male", "Male"}; %#ok<CLARRSTR>

% NOTE - For consistency purposes, let's use this dictionary to settle on
% colours for males and females. Later on (low priority) we can design a
% gui-function that prompts the user to select their own colour scheme.
ColorDict = dictionary(["Male", "Female"], ...
    {[0, 75, 75]./255, [245, 173, 82]./255});

AnimalReg = {IMOE.IMO_Data.Animals, IMOL.IMO_Data.Animals}; % Animal index per neuron
RF = 30; % FPS
EpochNames = {"Early", "Late"}; %#ok<CLARRSTR> 
GroupBy = Sexes;
Visualise = true;
Y_Label = "High movement events / second";

FreqEvents = F_Compare_HM_Events(Data1, Data2, LowBound, HighBound, ...
    RF, GroupBy, AnimalReg, EpochNames, Visualise, Y_Label, ColorDict,char(save_path));

image_EVF = strcat(path_stats_folder, "HM_Events\Sample - Event Frequency.pdf")
exportgraphics(gcf,image_EVF,"ContentType", "vector")
report_EVF = strcat(path_stats_folder, "HM_Events\Sample_HM_Freq.mat")
save(report_EVF, "FreqEvents")

Y_Label = "Mean event length";
[LengthEvents, Ls, Events] = F_Compare_HM_Lengths(Data1, Data2, LowBound, HighBound, ...
    RF, GroupBy, AnimalReg, EpochNames, Visualise, Y_Label, ColorDict,char(save_path));
Y_Label = "Mean event length";
image_EVL = strcat(path_stats_folder, "HM_Lenghts\Sample - Event Length.pdf")
exportgraphics(gcf,image_EVL,"ContentType", "vector")
report_EVL = strcat(path_stats_folder, "HM_Lenghts\Sample_HM_Len.mat")
save(report_EVL, "LengthEvents")
% 0.675

close(gcf)

%% Visualising coordination between motor and neuronal datasets
Animal = 2
flu = nansum(IMOL.IMO_Data.FiltTraces(IMOL.IMO_Data.AIX == Animal, :), 1);
flu = (flu - nanmedian(flu))/nanmedian(flu);

s1 = subplot(2, 1, 1);
plot(flu, 'Color', 'k')
ylabel("\DeltaF/F_0")
mov = MovL(IMOL.IMO_Data.Animals == Animal, :);

s2 = subplot(2, 1, 2);


plot(mov, 'Color', 'r')
ylabel('Movement score')
yticks([]);
yticklabels([]);
linkaxes([s1, s2], 'x')
