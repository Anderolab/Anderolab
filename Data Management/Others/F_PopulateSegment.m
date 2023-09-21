function [Peaks] = F_PopulateSegment(Segment, Factor)

peaks_ = randsample(Segment, floor(length(Segment)/Factor));
Peaks = zeros(size(Segment));
Peaks(peaks_) = 1;
area(Peaks)



end

