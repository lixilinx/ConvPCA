function B = PolyMatH(A)
% Hermitian of polynomial matrix
B = conj(permute(A, [2,1,3]));
B = B(:,:,end:-1:1);