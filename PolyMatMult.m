function C = PolyMatMult(A, B)
% polynomial matrices multiplication:
% (A0 + A1*z^-1 + A2*z^-2 + ...) * (B0 + B1*z^-1 + B2*z^-2 + ...)
% 
% to verify the implementation:
%   a = randn(1,1,3);
%   b = randn(1,1,5);
%   squeeze(PolyMatMult(a,b)) - conv(squeeze(a), squeeze(b)) % should be 0
%
[I,J,M] = size(A); [J,K,N] = size(B);
C = zeros(I,K,M+N-1);
for i = 0 : M+N-2
    for j = max(0, i+1-N) : min(M-1, i)
        C(:,:,i+1) = C(:,:,i+1) + A(:,:,j+1)*B(:,:,i-j+1);
    end
end