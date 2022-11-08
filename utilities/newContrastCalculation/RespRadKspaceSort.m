function spokePhases = RespRadKspaceSort(nFE,Times, RespWaveTimes,RespPhases)

currResp = 1;
count = 1;
spokePhases = zeros(nFE,1); % respiratory
lastRespphase = 1;
for i = 1:nFE
       RespTime = linspace(RespWaveTimes(currResp),RespWaveTimes(currResp+1),RespPhases+1);
       currTime = Times(i);
       if currTime > RespWaveTimes(currResp+1)
           currResp = currResp+1;
           if currResp > length(RespWaveTimes)
               currResp = length(RespWaveTimes);
           end
       end
       
       Respphase = discretize(currTime,RespTime); if isnan(Respphase); Respphase = 1; end
       spokePhases(i,1) = Respphase;
       if lastRespphase ~= Respphase
           count = count+1;
       end
       lastRespphase = Respphase;
end