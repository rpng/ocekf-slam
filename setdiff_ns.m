function [A] = setdiff_ns(A,B)
% Returns the values in A that are not in B, preserving the order of A.
% A and B are expected to be row vectors.

A = A(~ismembc(A,sort(B)));

if ~isempty(A)
    [T,R] = sort(A);
    I = [true diff(T)~=0];
    T = T(I);
    R = R(I);
    [I,I] = sort(R);
    A = T(I);
end 