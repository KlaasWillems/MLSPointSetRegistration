function lambda = KDE(X0, kernel_type, lambdas)
%KDE: this function implements kernel density estimation

epsilon = 1e-10; %I ensure that the density is sufficiently far from zero

%% initialization
% indx_struct=CG2.build_indx_struct(mesh);
% Mass = CG2.assembler.element_matrices(ones(size(struct_assembler.el.xq,1),1),struct_assembler.el,indx_struct.el,'mass');
N=size(X0,1);
if exist('lambdas','var')==0
    lambdas = N^(1/5)*logspace(-3,1,50);
end
Ncv=length(lambdas);

%%
switch kernel_type
    case 'Gaussian'
        phi = @(r)  exp(-r.^2);
        
    case 'matern'
        phi = @(r)  (1 + r + 1/3*r.^2) .*exp(-r);

    case 'multiquadric'
        phi= @(r) (1 + r.^2).^(-10);       
   
    otherwise	% default case
		error ('unknown kernel type')
end

if Ncv>1    
    CVs=zeros(Ncv,1);
    rcv = MeshfreePSR.util.stat.build_r_dr(X0,X0);
    %leave-one-out validation (this one tends to favor excessively large parameter values)
    for i=1:Ncv
        for j=1:N
            I = setdiff(1:N,j);
            rho = 1/lambdas(i)*mean(phi(rcv(j,I)/lambdas(i)),2); %I don't care about the normalization constants...
            CVs(i) = CVs(i) + 1/N*log(rho);
        end
    end
    [~,i]=max(CVs);
    lambda=lambdas(i);
else
    lambda=lambdas;
end

end