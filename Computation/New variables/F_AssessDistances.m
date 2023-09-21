function [Distance] = F_AssessDistances(Point, Centroid)
%F_ASSESSDISTANCES Summary of this function goes here
%   Detailed explanation goes here
Distance = sqrt(sum((Point-Centroid).^2, 2));
end

