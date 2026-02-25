function [lambdas,vectp] = compute_eigs(Hess)
%EIG_V (this only works for symmetric matrices! Hess=[Hess11,Hess12,Hess22])

%% initialization
N = size(Hess,1);
I = abs(Hess(:,2))>1e-10; %the indices in I correspond to the case in which the Hessian is not diagonal
J = not(I); %the indices in I correspond to the case in which the Hessian is diagonal
 
%%
trH  = Hess(I,1)+Hess(I,3);
discrH =sqrt((Hess(I,1)-Hess(I,3)).^2+4*Hess(I,2).^2); %this way of computing the discrimant is more stable
lambdas(I,1) = trH/2 - discrH/2;
lambdas(I,2) = trH/2 + discrH/2;

lambdas(J,1) = Hess(J,1);
lambdas(J,2) = Hess(J,3);

vectp = ones(N,4);
vectp(I,2)= -(Hess(I,1)-lambdas(I,1))./Hess(I,2);
vectp(I,4)= -(Hess(I,1)-lambdas(I,2))./Hess(I,2);
vectp(J,:)= repmat([1,0,0,1],sum(J),1);

%%
vectp(:,1:2)=vectp(:,1:2).*(1./sqrt(vectp(:,1).^2+vectp(:,2).^2));
vectp(:,3:4)=vectp(:,3:4).*(1./sqrt(vectp(:,3).^2+vectp(:,4).^2));

end

