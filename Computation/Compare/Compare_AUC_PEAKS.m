

%%
function Compare_AUC_PEAKS(AUC_Comparison, PEAKS_Comparison)

% Obtener la lista única de animales
animals = unique(AUC_Comparison.Animals);

% Crear un diccionario de colores para cada animal
colors = lines(length(animals));



% Inicializar la figura
figure
subplot(1, 2, 1) % Gráfico 1 para Early
hold on

% Bucle para recorrer la lista de animales y graficar sus puntos Early
for i = 1:length(animals)
    animal = animals(i);

    % Indices de Early para el animal actual
    early_idx = AUC_Comparison.Animals == animal;

    % Graficar puntos Early
    h=scatter(AUC_Comparison.Epochs(early_idx, 1), PEAKS_Comparison.Epochs_peak(early_idx, 1), 'Marker', 'o', 'MarkerFaceColor', colors(i, :), 'MarkerEdgeColor', 'k');
    h.UserData=animal; 

    % Realizar ajuste lineal y graficar la línea
%     coeffs = polyfit(AUC_Comparison.Epochs(early_idx, 1), PEAKS_Comparison.Epochs_peak(early_idx, 1), 1);
%     x_fit = min(AUC_Comparison.Epochs(early_idx, 1)):max(AUC_Comparison.Epochs(early_idx, 1));
%     y_fit = polyval(coeffs, x_fit);
%     plot(x_fit, y_fit, 'Color', colors(i, :));
end

% % Ajuste lineal para Early
early_x = AUC_Comparison.Epochs(:, 1);
early_y = PEAKS_Comparison.Epochs_peak(:, 1);
valid_idx = ~isnan(early_x) & ~isnan(early_y);
early_x_clean = early_x(valid_idx);
early_y_clean = early_y(valid_idx);
early_p = polyfit(early_x_clean, early_y_clean, 1);
early_line_y = polyval(early_p, early_x_clean);
plot(early_x_clean, early_line_y, 'r-'); % Ajuste lineal en rojo

early_y_fit = polyval(early_p, early_x_clean);

% Calcular el coeficiente de determinación (R²)
ssres = sum((early_y_clean - early_y_fit).^2);
sstot = sum((early_y_clean - mean(early_y_clean)).^2);
r_squared = 1 - ssres / sstot;

% Mostrar el coeficiente de determinación en el gráfico
text(max(early_x_clean), max(early_y_clean), ['R^2 = ' num2str(r_squared)], 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');

% Configurar límites y etiquetas de los ejes
xlabel('Peaks')
ylabel('AUC')
title('Early')

% Crear una leyenda personalizada para Early
legend_labels = strcat('Animal ', string(animals));
legend_colors = colors;
num_labels = length(legend_labels);
dummy_plots = zeros(1, num_labels);
for i = 1:num_labels
    dummy_plots(i) = plot(NaN, NaN, 'o', 'MarkerFaceColor', legend_colors(i, :), 'MarkerEdgeColor', 'k', 'DisplayName', legend_labels{i});
end
legend(dummy_plots, legend_labels);

% Mantener el gráfico
hold off

subplot(1, 2, 2) % Gráfico 2 para Late
hold on

% Bucle para recorrer la lista de animales y graficar sus puntos Late
for i = 1:length(animals)
    animal = animals(i);

    % Indices de Late para el animal actual
    late_idx = PEAKS_Comparison.Animals == animal;

    % Graficar puntos Late
    h=scatter(AUC_Comparison.Epochs(late_idx, 2), PEAKS_Comparison.Epochs_peak(late_idx, 2), 'Marker', 'o', 'MarkerFaceColor', colors(i, :), 'MarkerEdgeColor', 'k');
    h.UserData=animal;

     % Ajuste lineal de los puntos Early
    % Realizar ajuste lineal y graficar la línea
%     coeffs = polyfit(AUC_Comparison.Epochs(late_idx, 2), PEAKS_Comparison.Epochs_peak(late_idx, 2), 1);
%     x_fit = min(AUC_Comparison.Epochs(late_idx, 2)):max(AUC_Comparison.Epochs(late_idx, 2));
%     y_fit = polyval(coeffs, x_fit);
%     plot(x_fit, y_fit, 'Color', colors(i, :));

    % Ajuste lineal de los puntos Late
    

end

% % Ajuste lineal para Late
late_x = AUC_Comparison.Epochs(:, 2);
late_y = PEAKS_Comparison.Epochs_peak(:, 2);
valid_idx_late = ~isnan(late_x) & ~isnan(late_y);
late_x_clean = late_x(valid_idx_late);
late_y_clean = late_y(valid_idx_late);
late_p = polyfit(late_x_clean, late_y_clean, 1);
late_line_y = polyval(late_p, late_x_clean);
plot(late_x_clean, late_line_y, 'r-'); % Ajuste lineal en rojo
% late_p = polyfit(late_x, late_y, 1);
% late_line_y = polyval(late_p, late_x);
% plot(late_x, late_line_y, 'r-'); % Ajuste lineal en rojo

% Calcular valores ajustados
late_y_fit = polyval(late_p, late_x_clean);

% Calcular el coeficiente de determinación (R²)
ssres = sum((late_y_clean - late_y_fit).^2);
sstot = sum((late_y_clean - mean(late_y_clean)).^2);
r_squared = 1 - ssres / sstot;

% Mostrar el coeficiente de determinación en el gráfico
text(max(late_x_clean), max(late_y_clean), ['R^2 = ' num2str(r_squared)], 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');



% Configurar límites y etiquetas de los ejes
xlabel('Peaks')
ylabel('AUC')
title('Late')

% Crear una leyenda personalizada para Late
legend_labels = strcat('Animal ', string(animals));
legend_colors = colors;
num_labels = length(legend_labels);
dummy_plots = zeros(1, num_labels);
for i = 1:num_labels
    dummy_plots(i) = plot(NaN, NaN, 'o', 'MarkerFaceColor', legend_colors(i, :), 'MarkerEdgeColor', 'k', 'DisplayName', legend_labels{i});
end
legend(dummy_plots, legend_labels);

% Mantener el gráfico
hold off
% Habilitar el modo de cursor de datos y asignar la función de actualización personalizada
dcm_obj = datacursormode(gcf);
set(dcm_obj, 'UpdateFcn', @myupdatefcn);

% Función de actualización personalizada para mostrar el número del animal
function txt = myupdatefcn(~, event_obj)
pos = get(event_obj, 'Position');
animal_number = get(event_obj.Target, 'UserData');
txt = {['X: ', num2str(pos(1))], ['Y: ', num2str(pos(2))], ['Animal: ', num2str(animal_number)]};
end

end