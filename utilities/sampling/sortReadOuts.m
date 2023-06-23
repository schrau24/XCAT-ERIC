function sortedROs = sortReadOuts(nPhases, contrastTiming, nROs, readoutTiming)

sortedROs = zeros(nROs,1);
sortedROs(1) = 1;
for ms = 1:nPhases
    currTime = contrastTiming(ms);
    if ms<nPhases
        endTime = contrastTiming(ms+1);
    else
        endTime = max(readoutTiming);
    end
    temp_ind = find(readoutTiming>currTime & readoutTiming<=endTime);
    sortedROs(temp_ind) = ms;
end