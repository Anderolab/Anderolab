function [Glob_Mov, EdgePixels] = F_GetGlobalMovement(VideoObject)
% This fuction attemtps to extract global movement in the visual field,
% aming to generate a gross estimation of the movement of the animal in the
% visual field.

% Performance is low when animals are moving in a complex background, such
% as the home cage

%% Computing the edges and estimating movement
f_1 = read(VideoObject, 1);
last_Edge = F_GetEdge(rgb2gray(f_1(:, :, :)), .5);
Glob_Mov = 1:size(VideoObject.NumFrames, 4);
Glob_Mov(1) = 0;
EdgePixels(1) = sum(last_Edge==1, "all")/(size(last_Edge, 1)*size(last_Edge, 2));

% Looping through each frame
for fr = 2:VideoObject.NumFrames
    frame = read(VideoObject, fr);
    % 1 - Converting the frame to black and white and extracting the edges
    Edge = F_GetEdge(rgb2gray(frame(:, :, :)), .5);


    % 2 - Extracting and saving the edge differences
    edge_pixels = sum(Edge==1, "all");
    Glob_Mov(fr) = sum((Edge - last_Edge)==1, 'all')/edge_pixels;
    EdgePixels(fr) = ...
        edge_pixels/(size(Edge, 1)*size(Edge, 2));
    

    % Updating the previous edge
    last_Edge = Edge;


end
end

