function [y,x] = VanDerCorput(N)
% [y,x] = VanDerCorput(N);
%
% Van der Corput sequence for N points
% Generates a van der Corput sequence
dx = 1/N;
x = dx/2:dx:1-dx/2;
y = zeros(size(x));
for n=1:N
    nb = dec2bin(n);
    L = length(nb);
    nb = fliplr(nb);
    num = bin2dec(nb)/2^L;
    y(n) = num;
end