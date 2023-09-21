function [Table1,Table2,MdlLinear,coeffs]=Compare_Neurons(Data1,Data2,RF)
%%
dat1_prtl=prctile(Data1.RawTraces,5,2);
dat2_prtl=prctile(Data2.RawTraces,5,2);
incF1=(Data1.RawTraces-dat1_prtl)./(abs(dat1_prtl));
incF2=(Data2.RawTraces-dat2_prtl)./(abs(dat2_prtl));
data_filt={Data1.FiltTraces,Data2.FiltTraces};
Data = {incF1,incF2};
Data_sex={Data1.SEX,Data2.SEX};
Data_AIX={Data1.AIX,Data2.AIX};
numDataSets = length(Data);

for d_ix = 1:numDataSets
    select_data = Data{d_ix};
    size(select_data, 1);
    select_data_filt=data_filt{d_ix};
    PEAKS = zeros(size(select_data,1),1);
    AMPL = zeros(size(select_data,1),1);
    AUC = zeros(size(select_data,1),1);
    for row_ix = 1:size(select_data,1)

        peak_data = islocalmax(select_data_filt(row_ix,:), 2);
        PEAKS(row_ix) = (nansum(peak_data)/sum(~isnan(select_data(row_ix,:))))/RF; % cantidad de picos - Entre cantidad frames no nan
        AMPL(row_ix) = (nansum(peak_data .* select_data(row_ix,:))/nansum(peak_data)); % suma de amplitudes entre numero de picos
        AUC(row_ix) = (trapz(select_data(row_ix, ~isnan(select_data(row_ix, :))), 2)/sum(~isnan(select_data(row_ix,:))))/RF; % área bajo la curva entre longitud
    end
    SEX=Data_sex{d_ix}.';
    AIX=Data_AIX{d_ix}.';
    T = table(PEAKS,AMPL,AUC,SEX,AIX)
   
    size(T)
    if d_ix == 1
        Table1 = T
    else
        Table2 = T
    end
end
figure
hold on

scatter3(Table1.PEAKS(Table1.SEX == 'Male'), Table1.AMPL(Table1.SEX == 'Male'), Table1.AUC(Table1.SEX == 'Male'), 'k','filled')  % Machos en rojo
scatter3(Table1.PEAKS(Table1.SEX == 'Female'), Table1.AMPL(Table1.SEX == 'Female'), Table1.AUC(Table1.SEX == 'Female'), 'b','filled')  % Hembras en azul

% Graficar datos de Table2 en azul
% Distinguir entre machos y hembras en Data2
scatter3(Table2.PEAKS(Table2.SEX == 'Male'), Table2.AMPL(Table2.SEX == 'Male'), Table2.AUC(Table2.SEX == 'Male'), 'g','filled')  % Machos en verde
scatter3(Table2.PEAKS(Table2.SEX == 'Female'), Table2.AMPL(Table2.SEX == 'Female'), Table2.AUC(Table2.SEX== 'Female'), 'm','filled')  % Hembras en morado


% Añadir etiquetas a los ejes
xlabel('PEAKS')
ylabel('AMPL')
zlabel('AUC')

% Añadir una leyenda
legend('Male Early', 'Female Early', 'Male Late', 'Female Late')

hold off

% Concatena los datos de Table1 y Table2
AllData = [Table1; Table2];
AllData(:,4) = [];

% Define las etiquetas de las clases en función del estado y del sexo
labels = [strcat('Early_', Data1.SEX'); strcat('Late_', Data2.SEX')];
labels_categorical = categorical(labels);


% Realiza el análisis discriminante lineal
MdlLinear = fitcdiscr(AllData, labels_categorical);

% Muestra el modelo en consola
MdlLinear

% Obtener las puntuaciones discriminantes
[~,Scores] = MdlLinear.predict(AllData);
coeffs=MdlLinear.Coeffs;

% Crear un gráfico de caja y bigotes de las puntuaciones por grupo
figure
boxplot(Scores', labels_categorical, 'Notch', 'on', 'Labels', categories(labels_categorical));
xlabel('Group')
ylabel('Discriminant Score')
title('Boxplot of discriminant scores for each group')
% group1 = 'Early_Male';
% group2 = 'Early_Female';

% index = find((strcmp({MdlLinear.Coeffs.Class1}, group1) & strcmp({MdlLinear.Coeffs.Class2}, group2)) | ...
%              (strcmp({MdlLinear.Coeffs.Class1}, group2) & strcmp({MdlLinear.Coeffs.Class2}, group1)))
% % MdlLinear.Coeffs(index).Linear
% % constCoeffs = MdlLinear.Coeffs(index).Const
% Supongamos que Scores es un array cell, donde cada celda contiene las puntuaciones de un grupo
indices_Early_M = find(labels_categorical== 'Early_Male');
indices_Late_M = find(labels_categorical== 'Late_Male');
indices_Early_F = find(labels_categorical== 'Early_Female');
indices_Late_F = find(labels_categorical== 'Late_Female');
% Separar las puntuaciones por grupo
Scores_Early_M = Scores(indices_Early_M);
Scores_Late_M = Scores(indices_Late_M);
Scores_Early_F = Scores(indices_Early_F);
Scores_Late_F = Scores(indices_Late_F);

% Realizar el test de Kruskal-Wallis
[p, tbl, stats] = kruskalwallis([Scores_Early_M; Scores_Late_M; Scores_Early_F; Scores_Late_F], ...
    [repmat({'Early_Male'}, numel(Scores_Early_M), 1); repmat({'Late_Male'}, numel(Scores_Late_M), 1);...
    repmat({'Early_Female'}, numel(Scores_Early_F), 1); repmat({'Late_Female'}, numel(Scores_Late_F), 1)]);
disp(stats)
% Mostrar el valor p
disp(p)


% Realiza el test de Tukey
c = multcompare(stats,'CType','tukey-kramer');


