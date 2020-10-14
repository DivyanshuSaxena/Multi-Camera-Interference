% Find a set of binary vectors for primitive polynomials. By Jongho Lee 6/28/2018

% Input:
% stageN: number of stages

% Output:
% PrimPoly: a set of binary vectors representing primitive polynomials.
% Each row means a polynomial. e.g. [1 1 0 1] means x^3 + x^2 + 1

function PrimPoly = findPrimPoly(stageN)

pr = primpoly(stageN, 'all', 'nodisplay');
PrimPoly = de2bi(pr, 'left-msb');