function d_hat = estimateDepth_ACO(d, c, N, e_s, e_a, e_i, f_mod, T)


% Time-of-flight
tau = 2*d/c;                                        



%% Get correlations
C1 = T*(e_s + e_a + N*e_i + e_s/2.*cos(2*pi*f_mod*tau));
C2 = T*(e_s + e_a + N*e_i - e_s/2.*sin(2*pi*f_mod*tau));
C3 = T*(e_s + e_a + N*e_i - e_s/2.*cos(2*pi*f_mod*tau));
C4 = T*(e_s + e_a + N*e_i + e_s/2.*sin(2*pi*f_mod*tau));



%% Add Poisson noise
C1 = poissrnd(C1);
C2 = poissrnd(C2);
C3 = poissrnd(C3);
C4 = poissrnd(C4);



%% Decode
phase_hat = atan2((C4-C2) , (C1-C3));
phase_hat(phase_hat<0) = phase_hat(phase_hat<0) + 2*pi;
d_hat = c/(4*pi*f_mod)*phase_hat;