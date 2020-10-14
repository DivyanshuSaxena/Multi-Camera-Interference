function itfAmnt = estItfAmnt(N, binarySeq, start, ONIdx)

itfAmnt = zeros(1, size(ONIdx, 2));
for n = 1 : N
    
    ONIdxItf = find(binarySeq(n+1, :) == 1);
    startOneCam = start(n, 1);
    
    
    % Find interfered slots and interfering amount
    [itfIdx, itfOrderMsr, itfOrderItf] = intersect(ONIdx, ONIdxItf);
    itfAmnt(1, itfOrderMsr) = itfAmnt(1, itfOrderMsr) + 1 - abs(startOneCam);
    
    
    % Do one more time for other overlapping slots
    if startOneCam >= 0       
        ONIdxItf = ONIdxItf + 1;    
    else
        ONIdxItf = ONIdxItf - 1;
    end
    
    [itfIdx, itfOrderMsr, itfOrderItf] = intersect(ONIdx, ONIdxItf);
    itfAmnt(1, itfOrderMsr) = itfAmnt(1, itfOrderMsr) + abs(startOneCam);
      
end