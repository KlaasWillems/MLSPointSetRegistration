function coeff = labfm(dxji, dyji, hi, order)
% xji = [dx1i dx2i dx3i ]

assert(order >= 2); % We need all first-order and second-order derivatives.

p = (order*order + 3*order)/2;

Xt = MLS2d.LABFMUtil.getMonomials(dxji, dyji, order);
Wt = MLS2d.LABFMUtil.getBasisFunctions(dxji, dyji, hi, order);

M = Xt * Wt;

C = zeros(p, 5);
C(1:5, 1:5) = eye(5);

phi = M\C;

coeff = Wt*phi;

end