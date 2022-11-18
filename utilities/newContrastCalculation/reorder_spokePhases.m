function [GT_Resp] = reorder_spokePhases(spokePhases)
% ord1 = [1 2 3 4 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20];
ord2 = [1 1 1 2 3 4 5 6 7 8 9 8 7 6 5 4 3 2 1 1];

tmpSpokePhases = zeros(size(spokePhases));
for i = 1:20
    idx = find(spokePhases == i);
    tmpSpokePhases(idx) = ord2(i);
end
 
% flip and normalize
tmpSpokePhases = -tmpSpokePhases;
GT_Resp = tmpSpokePhases+9;