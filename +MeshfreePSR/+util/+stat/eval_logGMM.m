function logrho = eval_logGMM(xeval, GMM)
% Source: https://gregorygundersen.com/blog/2020/02/09/log-sum-exp/

% Optimized for large kappa (e.g., 4000), small dim = 2

[kappa, dim] = size(GMM.mu);
N = size(xeval,1);
tildexi = zeros(N, kappa);

X = xeval;          % (N×2)
MU = GMM.mu;        % (kappa×2)
omega = GMM.ComponentProportion(:);
const = (2*pi)^(dim/2);

% Precompute inverse covariances and constants
invSigma = zeros(2,2,kappa);
logw = zeros(kappa,1);

for i = 1:kappa
    S = GMM.Sigma(:,:,i);
    detS = S(1,1)*S(2,2) - S(1,2)*S(2,1);
    invSigma(:,:,i) = (1/detS) * [ S(2,2), -S(1,2); -S(2,1), S(1,1) ];
    logw(i) = log(omega(i) / (sqrt(detS) * const));
end

% Vectorized across N, scalar across kappa → fastest for big kappa
for i = 1:kappa
    dx = X - MU(i,:);                     % (N×2)
    Qi = sum((dx * invSigma(:,:,i)) .* dx, 2);
    tildexi(:,i) = -0.5 * Qi + logw(i);
end

% Log-sum-exp
c = max(tildexi, [], 2);
logrho = c + log(sum(exp(tildexi - c), 2));

end