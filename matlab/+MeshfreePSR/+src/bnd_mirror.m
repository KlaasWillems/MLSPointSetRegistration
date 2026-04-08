function [xnew, ynew, bnd_normal, d] = bnd_mirror(p,fd)
% Mirror a point across a boundary using an orthogonal projection

step = 1e-8;

pproj = bnd_proj(p, fd, step);
pnew = 2*pproj - p;

% Mirrored points
xnew = pnew(:,1);
ynew = pnew(:,2);

% Outward facing normal at boundary
bnd_normal = p - pproj;
bnd_normal = bnd_normal ./ vecnorm(bnd_normal, 2, 2);

% Distances to boundary
d = fd(p);
end





