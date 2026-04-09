from math import sqrt 
import torch
import scipy.io
import numpy as np

class SimulationParameters:
    def __init__(self, dt, num_steps, weightsXinf, mu_Xinf, inv_covariances_T, det_T, R):
        self.dt = dt
        self.num_steps = num_steps
        self.weightsXinf = weightsXinf
        self.mu_Xinf = mu_Xinf
        self.inv_covariances_T = inv_covariances_T
        self.det_T = det_T
        self.R = R

# ─── Function to compute log-density gradient of final GMM ─────
def compute_log_density_gradient(positions, weights, means, inv_covariances, det, device):
    diff = positions.unsqueeze(1) - means.unsqueeze(0)          # (N, K, 2)
    exponent = -0.5 * torch.einsum('nki,kij,nkj->nk', diff, inv_covariances, diff)
    
    log_norm = -0.5 * (2 * torch.log(torch.tensor(2 * torch.pi, device=device)) + torch.log(det))
    log_densities = log_norm.unsqueeze(0) + exponent
    
    # Log-sum-exp trick
    max_log = torch.max(log_densities, dim=1, keepdim=True)[0]
    log_p_final = max_log + torch.logsumexp(log_densities - max_log, dim=1, keepdim=True)
    
    densities = torch.exp(log_densities - max_log)
    p_final = torch.sum(weights.unsqueeze(0) * densities, dim=1, keepdim=True)
    
    grad_components = -torch.einsum('nki,kij->nkj', diff, inv_covariances)
    grad_components *= densities.unsqueeze(2)
    
    grad_p_final = torch.einsum('k,nkj->nj', weights, grad_components)
    grad_log_p = grad_p_final / (p_final + 1e-10)
    
    return grad_log_p

def loadMATLABGMM(filename, varName, device):

    # Load the .mat file
    mat = scipy.io.loadmat(filename, struct_as_record=False, squeeze_me=True)
    # print(mat.keys())
    gmm = mat[varName]

    # Extract parameters
    weights = np.array(gmm.ComponentProportion)   # shape: (K,)
    means = np.array(gmm.mu)                      # shape: (K, D)
    covariances = np.array(gmm.Sigma)             # shape: (D, D, K) or (K, D, D)

    # Ensure covariance shape is (K, D, D)
    if covariances.shape[0] != weights.shape[0]:
        covariances = np.transpose(covariances, (2, 0, 1))

    # Convert to PyTorch tensors
    weights_torch = torch.tensor(weights, dtype=torch.float32, device=device)
    means_torch = torch.tensor(means, dtype=torch.float32, device=device)
    covariances_torch = torch.tensor(covariances, dtype=torch.float32, device=device)

    # Print shapes
    # print("Weights:", weights_torch)        # (K,)
    # print("Means:", means_torch)            # (K, D)
    # print("Covariances:", covariances_torch)  # (K, D, D)
    return weights_torch, means_torch, covariances_torch

def sample_gmm_efficient(weights, means, covs, n_samples):
    K, D = means.shape

    cat = torch.distributions.Categorical(weights)
    comp_idx = cat.sample((n_samples,))  # (N,)

    samples = torch.zeros(n_samples, D, device=means.device)

    for k in range(K):
        mask = comp_idx == k
        num_k = mask.sum()
        if num_k > 0:
            mvn = torch.distributions.MultivariateNormal(
                means[k], covariance_matrix=covs[k]
            )
            samples[mask] = mvn.sample((num_k,))

    return samples

# ─── Main simulation loop ─────
def doSimulation(positions, masks, avg_positions, simParams, device):
    for i in range(simParams.num_steps):
        
        # Compute drift = ∇ log p_final(x)
        drift = compute_log_density_gradient(
            positions,
            simParams.weightsXinf,
            simParams.mu_Xinf,
            simParams.inv_covariances_T,
            simParams.det_T, 
            device
        )
        
        # Add noise
        noise = sqrt(2) * sqrt(simParams.dt) * torch.randn_like(positions)
        new_positions = positions + drift * simParams.dt + noise
        
        # Reflect boundary condition
        norms_new = torch.norm(new_positions, dim=1)
        inside = norms_new < simParams.R
        if inside.any():
            new_positions[inside] = 1.005 * simParams.R * new_positions[inside] / norms_new[inside].unsqueeze(1)
        
        positions = new_positions
        
        # Track average position of selected particles
        avg_pos = (masks.float().T @ positions) / masks.float().sum(dim=0).clamp(min=1).unsqueeze(1)
        avg_positions.append(avg_pos)
        
        if i % 50 == 0:  # reduced printing frequency
            drift_norms = torch.norm(drift, dim=1)
            print(f"t = {i * simParams.dt:.4f}  |  Max drift: {drift_norms.max().item():.4f}  |  Mean drift: {drift_norms.mean().item():.4f}")
    final_positions = positions.cpu().numpy()
    big_tensor = torch.stack(avg_positions)

    return final_positions, big_tensor