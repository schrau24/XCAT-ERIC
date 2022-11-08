function angular_segments = caspr_sampling_angular_segments(ang_map,rings)
% Make angular segments per concentric ring and return matrix

nSeg             = max(rings(:));
angular_segments = zeros(size(rings));
for n = 1 : nSeg
    idx = find(rings == n);
    tmp_angles = ang_map(idx);
    [~,idx_angles] = sort(tmp_angles);
    angular_segments(idx(idx_angles)) = 1 : numel(idx_angles);
end

% END
end