from math import sqrt 
import torch
import scipy.io
import numpy as np

class SimulationParameters:
    def __init__(self, dt, num_steps, weightsXinf, mu_Xinf, inv_covariances_T, det_T, R, small_region_radius):
        self.dt = dt
        self.num_steps = num_steps
        self.weightsXinf = weightsXinf
        self.mu_Xinf = mu_Xinf
        self.inv_covariances_T = inv_covariances_T
        self.det_T = det_T
        self.R = R
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

import torch

def line_circle_intersections(X1, X2, C, R):
    # X1, X2: (N, 2)
    # C: (2,) or (1,2)
    # R: float

    d = X2 - X1                    # (N, 2)
    f = X1 - C                     # (N, 2)

    a = (d * d).sum(dim=1)         # (N,)
    b = 2 * (f * d).sum(dim=1)     # (N,)
    c = (f * f).sum(dim=1) - R**2  # (N,)

    # Solve quadratic: t = (-b ± sqrt(b^2 - 4ac)) / (2a)
    disc = b**2 - 4*a*c            # (N,)
    sqrt_disc = torch.sqrt(disc)

    t1 = (-b - sqrt_disc) / (2*a)
    t2 = (-b + sqrt_disc) / (2*a)

    # Because we KNOW there is exactly one valid intersection
    # we choose the t that lies in (0, 1]
    t = torch.where((t1 > 0) & (t1 <= 1), t1, t2)

    # Compute intersection points
    P = X1 + d * t.unsqueeze(1)     # (N, 2)

    return P

def normAtIntersection(P, C):
    eps = 0.0
    dnOut = P - C
    dnNormOut = torch.norm(dnOut, dim=1, keepdim=True)
    return dnOut / (dnNormOut + eps);

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
    num_centers = centers.shape[0]
    avg_positions = []
    
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
        # norms_new = torch.norm(new_positions, dim=1)
        # inside = norms_new < simParams.R
        # if inside.any():
        #     new_positions[inside] = 1.005 * simParams.R * new_positions[inside] / norms_new[inside].unsqueeze(1)
        
        mask = computeMasks(positions, centers, simParams.small_region_radius)
        newMask = computeMasks(new_positions, centers, simParams.small_region_radius)
        particlesPerCircle = mask.sum(dim=0)
        pIn = newMask & (~mask)  # particles that entered a circle
        pOut = mask & (~newMask)  # particles that left a circle
        
        # Update position of circles
        for center in range(num_centers):
            if particlesPerCircle[center].item() > 0:
                # Old code
                # dOut = new_positions[pOut[:, center], :] - positions[pOut[:, center], :]
                # dIn = new_positions[pIn[:, center], :] - positions[pIn[:, center], :]
                # tempVel = dOut.sum(dim=0) - dIn.sum(dim=0)
                # centers[center, :] += tempVel / particlesPerCircle[center].item()

                intersectPointsOut = line_circle_intersections(positions[pOut[:, center], :], new_positions[pOut[:, center], :], centers[center, :], simParams.small_region_radius)
                intersectPointsIn = line_circle_intersections(positions[pIn[:, center], :], new_positions[pIn[:, center], :], centers[center, :], simParams.small_region_radius)
                nOut = normAtIntersection(intersectPointsOut, centers[center, :])
                nIn = normAtIntersection(intersectPointsIn, centers[center, :])
                dOut = new_positions[pOut[:, center], :] - positions[pOut[:, center], :]
                dIn = new_positions[pIn[:, center], :] - positions[pIn[:, center], :]
                tempVelOut = (dOut * nOut).sum(dim=1, keepdim=True)*nOut
                tempVelInt = (dIn * nIn).sum(dim=1, keepdim=True)*nIn
                centers[center, :] += (tempVelOut.sum(dim=0) - tempVelInt.sum(dim=0)) / particlesPerCircle[center].item()
            else:
                print(f"No particles in center: {center} at time step {i}")

        avg_positions.append(centers.clone())
        
        if i % 50 == 0:  # reduced printing frequency
            drift_norms = torch.norm(drift, dim=1)
            print(f"t = {i * simParams.dt:.4f}  |  Max drift: {drift_norms.max().item():.4f}  |  Mean drift: {drift_norms.mean().item():.4f}")
            
        # Update positions for next iteration
        positions = new_positions
        
    final_positions = positions.cpu().numpy()
    big_tensor = torch.stack(avg_positions)

    return final_positions, big_tensor