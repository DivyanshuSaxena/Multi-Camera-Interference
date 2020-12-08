function [itfAmntON, ONIdx] = estItfAmntAllCSMA(frac, binarySeq, start)

M = size(binarySeq, 2);

% For the primary camera: Choose the slots that would be ON tentatively.
primaryProb = binarySeq(1, :);
binarySeqItf = binarySeq(2:end, :);
itfAmnt = zeros(1, M);
ONIdx = ones(1, M);
count1 = 0; count2 = 0;

% The offsets of various interfering cameras in the range of [0+frac, 1+frac]
startSlot = zeros(size(start));
startSlot(start >= 0) = start(start >= 0);
startSlot(start < 0) = start(start < 0) + 1;

% Invariant: At the mth iteration, all the ON/OFF analysis for slots with
% start before the start of the mth primary camera slot has been done.
for m = 1 : M
    
    % In the mth iteration, we first check if the primary camera ON slot
    % can still be ON, by analyzing the mth slot of negON cameras and
    % (m-1)th slot of posON cameras.
    %
    % Further, we check the ON validity of mth slot of posON cameras and 
    % the (m+1)th slot of negON cameras.

    % Check if it will still remain ON
    posONprev = [];
    if m > 1
        posONprev = intersect(find(binarySeqItf(:, m-1) == 1), find(start >= 0));
    end
    negONcurr = intersect(find(binarySeqItf(:, m) == 1), find(start < 0));
    totalON = union(posONprev, negONcurr);
    
    % With a probability, sense else defer to the next slot.
    if ~isempty(totalON)
        % If the offset is within some threshold, still switch the primary
        % camera ON.
        [offset, ~] = min(startSlot(totalON));
        if offset < 0.95
            count1 = count1 + 1;
            ONIdx(m) = 0;
        else
            if primaryProb(m) == 1
                itfAmnt(1, m) = itfAmnt(1, m) + offset - frac;
            else
                count1 = count1 + 1;
                ONIdx(m) = 0;
            end
        end
    else
        if primaryProb(m) == 0
            count2 = count2 +1;
            ONIdx(m) = 0;
        end
    end
    
    posONcurr = intersect(find(binarySeqItf(:, m) == 1), find(start >= 0));
    negONnext = [];
    if m < M
        negONnext = intersect(find(binarySeqItf(:, m+1) == 1), find(start < 0));
    end
    tentativeON = union(posONcurr, negONnext);
    
    if ONIdx(m) == 1
        % All tentativeON shall have to be OFF
        binarySeqItf(posONcurr, m) = 0;
        binarySeqItf(negONnext, m+1) = 0;
    else
        % Maximum one of the tentativeON shall be ON - decided based on a
        % probability p of sensing to be actually done.
        if ~isempty(tentativeON)
            [~, ONIndex] = min(startSlot(tentativeON));
            
            % For all others in tentativeON, set them as OFF
            posONcurr(posONcurr == ONIndex) = [];
            if ~isempty(negONnext)
                negONnext(negONnext == ONIndex) = [];
            end
            
            binarySeqItf(posONcurr, m) = 0;
            binarySeqItf(negONnext, m+1) = 0;
        end
    end
end

count1, count2
ONIdx = find(ONIdx == 1);
itfAmntON = itfAmnt(ONIdx);