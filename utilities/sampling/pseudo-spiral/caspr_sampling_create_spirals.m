function caspr = caspr_sampling_create_spirals(rings,angular_segments)
%Combine the rings and angular segments 

Sl          = max(rings(:));
Ns          = max(angular_segments(:));
caspr       = zeros(Sl,Ns);

for s = 1 : Ns
    ang_idx = 1 + mod(s - 1 + round(linspace(1,Ns,Sl)),Ns);

    for p = 1 : Sl
        tmp_idx = find(rings == p);
        tmp_idx2 = find(angular_segments(tmp_idx) == ang_idx(p));
        if ~isempty(tmp_idx2)
            caspr(p,s) = tmp_idx(tmp_idx2);
        else
            caspr(p,s) = caspr(p-1,s);
        end
    end
end

% END
end