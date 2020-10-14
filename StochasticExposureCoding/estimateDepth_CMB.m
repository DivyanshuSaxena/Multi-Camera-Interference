function d_hat = estimateDepth_CMB(d, c, p, N, M, A, e_s, e_a, e_i, f_mod, T)



%% Parameters
tau = 2*d/c;                                        % time-of-flight
sampleN = size(e_s, 1);                             % number of samples
Tslot = T/M;                                        % slot integration time



%% Generate binary sequences
binarySeq = binornd(1, p, N+1, M);                  % N+1 by M matrix: 1st row for the primary camera, 2nd~Nth rows for the interfering cameras



%% Generate random starting time for the interfering cameras
start = 2*rand(N, 1) - 1;                           % -1.0 ~ 1.0



%% Find the ON slots of the primary camera
ONIdx = find(binarySeq(1, :) == 1);
M_ON = size(ONIdx, 2);                              % number of the ON slots



%% Estimate interference amount due to the interfering cameras
itfAmnt = estItfAmnt(N, binarySeq, start, ONIdx);   % 1 by M_ON vector where each element contains the interference amount



%% Get correlation values
C1 = zeros(sampleN, M_ON);
C2 = zeros(sampleN, M_ON);
C3 = zeros(sampleN, M_ON);
C4 = zeros(sampleN, M_ON);

for m = 1 : M_ON
    
    C1(:, m) = Tslot*(A*e_s + e_a + itfAmnt(1, m)*A*e_i + A*e_s/2.*cos(2*pi*f_mod.*tau));
    C2(:, m) = Tslot*(A*e_s + e_a + itfAmnt(1, m)*A*e_i - A*e_s/2.*sin(2*pi*f_mod.*tau));
    C3(:, m) = Tslot*(A*e_s + e_a + itfAmnt(1, m)*A*e_i - A*e_s/2.*cos(2*pi*f_mod.*tau));
    C4(:, m) = Tslot*(A*e_s + e_a + itfAmnt(1, m)*A*e_i + A*e_s/2.*sin(2*pi*f_mod.*tau));
end



%% Add Poisson noise
C1 = poissrnd(C1);
C2 = poissrnd(C2);
C3 = poissrnd(C3);
C4 = poissrnd(C4);



%% Sum correlation values
C1 = sum(C1, 2);
C2 = sum(C2, 2);
C3 = sum(C3, 2);
C4 = sum(C4, 2);



%% Decode
phase_hat = atan2((C4-C2) , (C1-C3));
phase_hat(phase_hat<0) = phase_hat(phase_hat<0) + 2*pi;
d_hat = c/(4*pi*f_mod)*phase_hat;



