function [] = F_ViewDistribution(Variable, xlab, ylab)
%F_VIEWDISTRIBUTION Summary of this function goes here
%   Detailed explanation goes here
subplot(1, 3, 3)
histogram(Variable, "FaceAlpha", 1, "FaceColor", 'k');
ylabel(xlab)
xl = xlim();
set(gca,'view',[90 -90])
xticks([]);
subplot(1, 3, [2, 1])
boxplot(Variable)
ylim(xl)
xticks([])
ylabel(ylab)

end

