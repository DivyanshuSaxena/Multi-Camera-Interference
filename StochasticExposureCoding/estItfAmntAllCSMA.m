function [ONIdx] = estItfAmntAllCSMA(binarySeq, start)

M = size(binarySeq, 2);

% For the primary camera: Choose the slots that would be ON tentatively.
ONIdx = binarySeq(1, :);
binarySeqItf = binarySeq(2:end, :);

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

    if ONIdx(m) == 1
        % Check if it will still remain ON
        posONprev = [];
        if m > 1
            posONprev = intersect(find(binarySeqItf(:, m-1) == 1), find(start >= 0));
        end
        negONcurr = intersect(find(binarySeqItf(:, m) == 1), find(start < 0));
        totalON = union(posONprev, negONcurr);
        if ~isempty(totalON)
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
        if ~isempty(tentativeON)
            % Only one of the tentativeON shall be ON
            ONOffset = min(startSlot(tentativeON));
            ONIndex = find(startSlot == ONOffset);

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

ONIdx = find(ONIdx == 1);