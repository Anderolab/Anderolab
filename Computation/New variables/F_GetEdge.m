function [BWs] = F_GetEdge(Frame, FudgeFactor)

% Estimating the edge threshold via SOBEL method
% (SOBEL due to speed)
[~, tshld] = edge(Frame, 'sobel');

% And calculating the edges
BWs = edge(Frame, 'sobel', tshld*FudgeFactor);
end

