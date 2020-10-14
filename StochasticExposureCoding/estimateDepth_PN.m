function d_hat = estimateDepth_PN(d, c, N, ns, na, ni, f_mod, T, stageN, sampleNperBit)


%% Parameters
tau = 2*d/c;
T_mod = 1/f_mod;
bitN = 2^stageN - 1;
a = 2*bitN/(bitN + 1);
tDummy = 0 : T_mod/(sampleNperBit*bitN - 1) : T_mod;



%% Generate N+1 m-sequences
% Generate primitive polynomials. The number of primitive polynomials
% should be larger than N+1
PrimPoly = findPrimPoly(stageN);

if (size(PrimPoly, 1) < N+1)
    error('The number of primitive polynomials is less than N+1');
end


% Generate m-sequences with primitive polynomials
mSeq = zeros(N+1,  bitN);

initCond = zeros(1, stageN);
initCond(1, 1) = 1;
for n = 1 : N+1
   
    pnSeq = comm.PNSequence('Polynomial', PrimPoly(n, :), 'SamplesPerFrame', bitN, 'InitialConditions', initCond);
    mSeq(n, :) = pnSeq();
end



%% Generate a modulation function
modSet = repelem(mSeq, 1, sampleNperBit);
demod = modSet(1, :);                                       % Demodulation function

% Draw
% figure; plot(mod(1,:), 'r', 'lineWidth', 4)

% Autocorrelaiton
% autoCorr = computeCorrNumeric(mod(1, :), mod(1, :), T_mod/bitN/sampleNperBit);
% figure; plot(autoCorr, 'r', 'lineWidth', 4)



%% Delay
tauSample = round(tau/T_mod * size(modSet, 2));

% Measuring modulation function
modSet(1, :) = circshift(modSet(1, :), tauSample);

% Interfering modulation functions
for n = 1 : N
    
    delay = randi(size(modSet, 2) - 1);
    modSet(n+1, :) = circshift(modSet(n+1, :), delay);
end



%% Get correlations
modSet(1, :) = ns*modSet(1, :);
for n = 1 : N

    modSet(n+1, :) = ni*modSet(n+1, :);
end

mod = sum(modSet, 1);

C1 = 2*a*T/T_mod*trapz(tDummy, mod.*demod) + na*T*(bitN + 1)/bitN;

modSet(1, :) = circshift(modSet(1, :), -sampleNperBit);
mod = sum(modSet, 1);
C2 = 2*a*T/T_mod*trapz(tDummy, mod.*demod) + na*T*(bitN + 1)/bitN;

modSet(1, :) = circshift(modSet(1, :), -sampleNperBit);
mod = sum(modSet, 1);
C3 = 2*a*T/T_mod*trapz(tDummy, mod.*demod) + na*T*(bitN + 1)/bitN;

modSet(1, :) = circshift(modSet(1, :), -sampleNperBit);
mod = sum(modSet, 1);
C4 = 2*a*T/T_mod*trapz(tDummy, mod.*demod) + na*T*(bitN + 1)/bitN;



%% Add Poisson noise
C1 = poissrnd(C1);
C2 = poissrnd(C2);
C3 = poissrnd(C3);
C4 = poissrnd(C4);



%% Decode
tau_hat = T_mod/bitN*(C2 - C4)/(C1 + C2 - C3 - C4);
d_hat = tau_hat*c/2;



