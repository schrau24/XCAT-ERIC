% RMSSD         : RMS beat-to-beat differences (in ms) (6 is good)
% Constraint    : Multiplicative factor to bound HRV (0.14 is good)
% TraceLength   : Time of measurement (in seconds)
% BaseRR        : Base RR length

% RR            : Simulated RR intervals (CTG Trace)


function [RR] = generate_RR(RMSSD, Constraint, TraceLength, BaseRR)

clear HR
clear RR
RR(1) = BaseRR + RMSSD*randn(1,1);
while sum(RR) < TraceLength
    RR(end+1) = RR(end) + RMSSD*randn(1,1) - Constraint*(RR(end)-BaseRR);
end

% figure(2)
% plot(cumsum(RR)/1000, 60000./RR)

end