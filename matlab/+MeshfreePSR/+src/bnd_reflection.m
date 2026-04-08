function [xnew,ynew,pproj] = bnd_reflection(p,v,fd)

%% initialization
N=size(p,1);
tol=1e-8;

%%
pguess=p+v;
d=fd(pguess);
I = d<0;
% sprintf('Projecting %i points', sum(~I))
pnew=zeros(N,2);
pnew(I,:)=pguess(I,:); %if the new points are inside the domain we do not do anything
pproj = bnd_proj(pguess(~I,:),fd,tol);%otherwise, we compute the projection onto the boundary
pnew(~I,:)= 2*pproj - pguess(~I,:);

xnew=pnew(:,1);
ynew=pnew(:,2);
end





