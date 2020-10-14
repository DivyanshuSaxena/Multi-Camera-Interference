function noClshIdx_hat = checkClash(C1, C2, C3, C4, ONIdx, k)

C1sum = sum(C1, 1);
C2sum = sum(C2, 1);
C3sum = sum(C3, 1);
C4sum = sum(C4, 1);

o = (C1sum + C2sum + C3sum + C4sum);             % 1 x M_ON vector
o_min = min(o);
o_mean = o_min + k^2/2 + sqrt(k^2*o_min + k^4/4);
o_clsh = o_mean + k*sqrt(o_mean);

noClshIdx_hat = ONIdx(o < o_clsh);