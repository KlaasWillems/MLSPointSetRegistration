function [r, drdx, d2rdx2] = build_r_dr(X1,X2)

%% step 0: initialization
N1 = size(X1,1);
N2 = size(X2,1);
dim = size(X2,2);

%% step 1: I build r
norms1 = sum(X1.^2,2);
norms2 = sum(X2.^2,2);

mat1 = repmat(norms1,1,N2);
mat2 = repmat(norms2',N1,1);

r = real( sqrt(mat1 + mat2 - 2*X1*X2') );	% full distance matrix

%% step 2: I build drdx and drdy
if nargout>1
    drdx = zeros(N1,N2,dim);
    for d=1:dim
        drdx(:,:,d) = (1./(r+1e-16)).*(repmat(X1(:,d),1,N2) - repmat(X2(:,d)',N1,1));
    end
end
if nargout > 2
    % drdxdy
    d2rdx2 = zeros(N1, N2, dim+1);
    for d=1:dim
        d2rdx2(:,:,d) = (1./((r.^3)+1e-16)).*(repmat(X1(:,d),1,N2) - repmat(X2(:,d)',N1,1)).^2;
    end
    d2rdx2(:,:,dim+1) = (-1./((r.^3)+1e-16)).*(repmat(X1(:,1),1,N2) - repmat(X2(:,1)',N1,1)).*(repmat(X1(:,2),1,N2) - repmat(X2(:,2)',N1,1));
end
end