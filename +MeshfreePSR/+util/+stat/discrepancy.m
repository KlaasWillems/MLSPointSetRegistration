function d = discrepancy(x, y)
% Given x, y the x and y coordinates of points in the unit square, compute
% the discrepancy.

N = length(x);
ds = zeros(N, 1);
eps = 1e-8;

for i = 1:N
    xm = x(i)-eps;
    ym = y(i)-eps;
    xp = x(i)+eps;
    yp = y(i)+eps;

    d1 = xm*ym - sum( x <= xm & y <= ym)/N;
    d2 = xp*yp - sum( x <= xp & y <= yp)/N;
    ds(i) = max(d1, d2);
end

d = max(ds);
end