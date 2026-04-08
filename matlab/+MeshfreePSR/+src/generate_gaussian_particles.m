function [xb,yb,rhof] = generate_gaussian_particles(N,mu,Sigma)

%% initialization
[V, D] = eig(Sigma);
Vstar = sqrt(D)*V';

%% generate the initial condition for the standard Gaussian
[Y, X] = VanDerCorput(N);
x = erfinv(2*X(:)-1)*sqrt(2);
y = erfinv(2*Y(:)-1)*sqrt(2);

%% transformation
xb=mu(1)+Vstar(1,1)*x+Vstar(2,1)*y;
yb=mu(2)+Vstar(1,2)*x+Vstar(2,2)*y;
rhof=mvnpdf([xb,yb], mu, Sigma);

end
