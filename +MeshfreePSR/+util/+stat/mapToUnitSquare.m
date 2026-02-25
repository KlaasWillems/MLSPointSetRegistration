function [xnew, ynew] = mapToUnitSquare(x, y, w, mux, muy, sigma)
% Maps the points x, y to the unit square
% The mapping of points is the map that maps f in the Lebesgue measure in
% [0, 1]^2. In this case, f is given by a weighted (with wi) sum of
% Gaussians with means mux ans muy and standard deviations sigmai.
Ng = length(w);
N = length(x);
xnew = zeros(N, 1);
ynew = zeros(N, 1);

gi = @(t, si) exp(-(t^2)/(2*si*si))/(sqrt(2*pi)*si);
phi = @(t, si) (1 + erf( t/(sqrt(2)*si) ))/2;

for i = 1:N
    xi = x(i);
    yi = y(i);

    % Evaluate T1 and T2
    t1num = 0.0;
    t1denum = 0.0;
    t2 = 0.0;
    for j = 1:Ng
        val = gi(yi - muy(j), sigma(j))*w(j);
        t1num = t1num + val*phi(xi - mux(j), sigma(j));
        t1denum = t1denum + val;
        t2 = t2 + phi(yi - muy(j), sigma(j))*w(j);
    end
    xnew(i) = t1num/t1denum;
    ynew(i) = t2;
end

end