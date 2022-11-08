function caspr = caspr_sampling_shuffle_ga(caspr_uni,ang_map,Ns)
%Shuffle the caspr interleaves with the golden angle

if nargin < 3
    Ns = size(caspr_uni,2);
end

if size(caspr_uni,1) > 1
    theta_uni = pi + ang_map(caspr_uni(2,:));
else
    caspr = caspr_uni;
    return;
end

ga = (3. - sqrt(5.)) * pi;    % Golden angle in radians

% % tiny golden angles
% tau = (1+sqrt(5))/2;
% ga = pi / (tau + 7 - 1);

idx = [];
for n = 1 : Ns
    [~,idx(n)] = min(abs(mod(n*ga,2*pi)-theta_uni));
    theta_uni(idx(n)) = 1E06;
end

caspr = caspr_uni(:,idx);

% END
end