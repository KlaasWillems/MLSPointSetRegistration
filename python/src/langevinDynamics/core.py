from math import sqrt 
import torch
import scipy.io
import numpy as np

class SimulationParameters:
    def __init__(self, dt, num_steps, weightsXinf, mu_Xinf, inv_covariances_T, det_T, R, reflectiveObject, small_region_radius):
        self.dt = dt
        self.num_steps = num_steps
        self.weightsXinf = weightsXinf
        self.mu_Xinf = mu_Xinf
        self.inv_covariances_T = inv_covariances_T
        self.det_T = det_T
        self.R = R
        self.reflectiveObject = reflectiveObject
        self.small_region_radius = small_region_radius

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

def computeMasks(positions, centers, radius):
    num_particles = positions.shape[0]
    num_centers = centers.shape[0]
    masks = torch.empty((num_particles, num_centers), dtype=torch.bool, device=positions.device)
    
    for j in range(num_centers):
        center = centers[j]
        dist_squared = torch.sum((positions - center) ** 2, dim=1)
        masks[:, j] = dist_squared < radius ** 2
    
    return masks

# ─── Main simulation loop ─────
def doSimulation(positions, centers, simParams, device):
    reflectionVar = 1.002
    num_centers = centers.shape[0]
    avg_positions = torch.empty((simParams.num_steps, num_centers, 2), dtype=torch.float32, device=device)
    amountOfParticlesPerCircle = torch.empty((simParams.num_steps, num_centers), dtype=torch.int32, device=device)
    
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
        
        # Reflect particles
        norms_new = torch.norm(new_positions, dim=1)
        inside = norms_new < simParams.R
        if inside.any() and simParams.reflectiveObject:
            new_positions[inside] = reflectionVar * simParams.R * new_positions[inside] / norms_new[inside].unsqueeze(1)
        
        mask = computeMasks(positions, centers, simParams.small_region_radius)
        driftOfCircles = compute_log_density_gradient(centers, simParams.weightsXinf, simParams.mu_Xinf, simParams.inv_covariances_T, simParams.det_T, device)  # (2,)
        particlesPerCircle = mask.sum(dim=0)
        
        # Update position of circles
        for center in range(num_centers):
            if particlesPerCircle[center].item() > 0:
                shifts = positions[mask[:, center], :] - centers[center, :]  # (num_particles_in_circle, 2)
                meanShiftVector = shifts.mean(dim=0)  # (, 2)
                centers[center, :] += simParams.dt*( driftOfCircles[center, :] - 4*meanShiftVector/(simParams.small_region_radius**2))
            else:
                print(f"No particles in center: {center} at time step {i}")
                
        # Reflect circles
        norms_new = torch.norm(centers, dim=1)
        inside = norms_new < simParams.R
        if inside.any() and simParams.reflectiveObject:
            centers[inside] = reflectionVar * simParams.R * centers[inside] / norms_new[inside].unsqueeze(1)

        amountOfParticlesPerCircle[i, :] = particlesPerCircle
        avg_positions[i, :, :] = centers.clone()

        if i % 50 == 0:  # reduced printing frequency
            drift_norms = torch.norm(drift, dim=1)
            print(f"t = {i * simParams.dt:.4f}  |  Max drift: {drift_norms.max().item():.4f}  |  Mean drift: {drift_norms.mean().item():.4f}")
            
        # Update positions for next iteration
        positions = new_positions
        
    final_positions = positions.cpu().numpy()

    return final_positions, avg_positions, amountOfParticlesPerCircle.cpu().numpy()

# ─── Main simulation loop ─────
def doSimulationTry(positions, centers, simParams, device):
    reflectionVar = 1.002
    num_centers = centers.shape[0]
    avg_positions = torch.empty((simParams.num_steps, num_centers, 2), dtype=torch.float32, device=device)
    amountOfParticlesPerCircle = torch.empty((simParams.num_steps, num_centers), dtype=torch.int32, device=device)
    
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
        new_positionsTilde = positions + noise
        
        # Reflect particles
        # norms_new = torch.norm(new_positions, dim=1)
        # inside = norms_new < simParams.R
        # if inside.any() and simParams.reflectiveObject:
        #     new_positions[inside] = reflectionVar * simParams.R * new_positions[inside] / norms_new[inside].unsqueeze(1)
        
        mask = computeMasks(new_positionsTilde, centers, simParams.small_region_radius)
        driftOfCircles = compute_log_density_gradient(centers, simParams.weightsXinf, simParams.mu_Xinf, simParams.inv_covariances_T, simParams.det_T, device)  # (2,)
        particlesPerCircle = mask.sum(dim=0)
        
        # Update position of circles
        for center in range(num_centers):
            if particlesPerCircle[center].item() > 0:
                meanBackwardDrift = (new_positionsTilde[mask[:, center], :] - positions[mask[:, center], :]).mean(dim=0) / simParams.dt
                centers[center, :] += simParams.dt*( driftOfCircles[center, :] + meanBackwardDrift/2)
            else:
                print(f"No particles in center: {center} at time step {i}")
                
        # Add drift
        new_positions = new_positionsTilde + drift * simParams.dt
                
        # Reflect circles
        # norms_new = torch.norm(centers, dim=1)
        # inside = norms_new < simParams.R
        # if inside.any() and simParams.reflectiveObject:
        #     centers[inside] = reflectionVar * simParams.R * centers[inside] / norms_new[inside].unsqueeze(1)

        amountOfParticlesPerCircle[i, :] = particlesPerCircle
        avg_positions[i, :, :] = centers.clone()

        if i % 50 == 0:  # reduced printing frequency
            drift_norms = torch.norm(drift, dim=1)
            print(f"t = {i * simParams.dt:.4f}  |  Max drift: {drift_norms.max().item():.4f}  |  Mean drift: {drift_norms.mean().item():.4f}")
            
        # Update positions for next iteration
        positions = new_positions
        
    final_positions = positions.cpu().numpy()

    return final_positions, avg_positions, amountOfParticlesPerCircle.cpu().numpy()