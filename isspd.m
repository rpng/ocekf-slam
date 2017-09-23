
function [t,R] = isspd(Sigma)
%ISPDS Test if a matrix is positive definite symmetric
% T = ISPDS(SIGMA) returns a logical indicating whether the matrix SIGMA is
% square, symmetric, and positive definite, i.e., it is a valid full rank
% covariance matrix.
%
% [T,R] = ISPDS(SIGMA) returns the cholesky factor of SIGMA in R.  If SIGMA
% is not square symmetric, ISPDS returns [] in R.

% Test for square, symmetric
%
[nSamples,m] = size(Sigma);
if (nSamples == m) & all(all(abs(Sigma - Sigma') < 10*eps(max(abs(diag(Sigma))))))
    % Test for positive definiteness
    [R,p] = chol(Sigma);
    if p == 0
        t = true;
    else
        t = false;
    end
else
    R = [];
    t = false;
end