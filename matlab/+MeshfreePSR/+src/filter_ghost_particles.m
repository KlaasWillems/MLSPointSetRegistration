function [xmirror, ymirror] = filter_ghost_particles(x, y, xmirror, ymirror, inputs_method)

% Ghost particles that lie inside the domain are removed by setting their
% coordinate to (Inf, Inf) (they won't be added as neighbours then).
dist = inputs_method.fd([xmirror, ymirror]);
xmirror(dist < 0) = Inf;
ymirror(dist < 0) = Inf;

% Ghost particles that lie farther than r away from their associated
% particle are also removed.
dist = sqrt((x - xmirror).^2 + (y - ymirror).^2);
xmirror(dist > inputs_method.Cr*2) = Inf;
ymirror(dist > inputs_method.Cr*2) = Inf;

end