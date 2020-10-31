function [itfAmntON, ONIdx] = estItfAmntCSMA(N, binarySeq, start)

m = size(binarySeq, 2);
itfAmnt = zeros(1, m);
ONIdx = 1:m;

for n = 1 : N
    
    ONIdxItf = find(binarySeq(n, :) == 1);
    startOneCam = start(n, 1);  % Offset > 0 implies interfering camera started after the camera
    
    if startOneCam < 0
        % Interfering camera already started.
        % Put an OFF for all interfering slots - the light won't be
        % transmitted in that scenario
        [itfIdx, itfOrderMsr, itfOrderItf] = intersect(ONIdx, ONIdxItf);
        ONIdx(itfOrderMsr) = [];
        
        % For the other overlapping slot, calculate the interference amount
        ONIdxItfNext = ONIdxItf - 1;
        [itfIdx, itfOrderMsr, itfOrderItf] = intersect(ONIdx, ONIdxItfNext);
        itfAmnt(1, itfOrderMsr) = itfAmnt(1, itfOrderMsr) + abs(startOneCam);
    end
    
    if startOneCam > 0
        % Interfering camera started after the primary camera
        % Put an OFF for all slots where the previous slot of itf camera
        % overlaps with a slot of primary camera
        ONIdxItfPrev = ONIdxItf + 1;
        [itfIdx, itfOrderMsr, itfOrderItf] = intersect(ONIdx, ONIdxItfPrev);
        ONIdx(itfOrderMsr) = [];
        
        % For the other overlapping slot, calculate the interference amount
        [itfIdx, itfOrderMsr, itfOrderItf] = intersect(ONIdx, ONIdxItf);
        itfAmnt(1, itfOrderMsr) = itfAmnt(1, itfOrderMsr) + 1 - abs(startOneCam);    
    end
      
end

itfAmntON = itfAmnt(ONIdx);