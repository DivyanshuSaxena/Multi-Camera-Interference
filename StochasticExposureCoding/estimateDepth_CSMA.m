function [d_hat, M_ON] = estimateDepth_CSMA(d, c, p, N, M, A, A_SEC, e_s, e_a, e_i, f_mod, T, frac)



%% Parameters
tau = 2*d/c;                                        % time-of-flight
sampleN = size(e_s, 1);                             % number of samples
Tslot = frac*T/M;                                   % slot integration time



%% Generate binary ON/OFF sequences
binarySeq = binornd(1, p, N, M);                  % N by M matrix: N rows for the interfering cameras


%% Generate random starting time for the interfering cameras
start = 2*rand(N, 1) - 1;                           % -1.0 ~ 1.0



%% Estimate interference amount due to the interfering cameras
[itfAmnt, ONIdx] = estItfAmntCSMA(N, binarySeq, start, frac);   % 1 by M vector where each element contains the interference amount
M_ON = size(itfAmnt, 2);
assert(size(itfAmnt, 2) == size(ONIdx, 2));



%% Get correlation values
C1 = zeros(sampleN, M_ON);
C2 = zeros(sampleN, M_ON);
C3 = zeros(sampleN, M_ON);
C4 = zeros(sampleN, M_ON);

for m = 1 : M_ON
    
    C1(:, m) = Tslot*(A*e_s + e_a + itfAmnt(1, m)*A_SEC*e_i + A*e_s/2.*cos(2*pi*f_mod.*tau));
    C2(:, m) = Tslot*(A*e_s + e_a + itfAmnt(1, m)*A_SEC*e_i - A*e_s/2.*sin(2*pi*f_mod.*tau));
    C3(:, m) = Tslot*(A*e_s + e_a + itfAmnt(1, m)*A_SEC*e_i - A*e_s/2.*cos(2*pi*f_mod.*tau));
    C4(:, m) = Tslot*(A*e_s + e_a + itfAmnt(1, m)*A_SEC*e_i + A*e_s/2.*sin(2*pi*f_mod.*tau));
end



%% Add Poisson noise
C1 = poissrnd(C1);
C2 = poissrnd(C2);
C3 = poissrnd(C3);
C4 = poissrnd(C4);



%% Check clash
% noClshIdx = noClshIdxGT;                          % Use ground truth

k = 2;
noClshIdx = checkClash(C1, C2, C3, C4, ONIdx, k);   % Use clash check algorithm
% [size(ONIdx,2), size(noClshIdx, 2)]



%% Extract non-clashed slots
[yesno, memberIdx] = ismember(noClshIdx, ONIdx);
C1 = C1(:, memberIdx);
C2 = C2(:, memberIdx);
C3 = C3(:, memberIdx);
C4 = C4(:, memberIdx);



%% Sum correlation values
C1 = sum(C1, 2);
C2 = sum(C2, 2);
C3 = sum(C3, 2);
C4 = sum(C4, 2);



%% Decode
phase_hat = atan2((C4-C2) , (C1-C3));
phase_hat(phase_hat<0) = phase_hat(phase_hat<0) + 2*pi;
d_hat = c/(4*pi*f_mod)*phase_hat;



