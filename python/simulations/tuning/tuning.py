import argparse
from pathlib import Path

import numpy as np
import scipy.io


DEFAULT_RADII = np.logspace(-2, 0, 10, dtype=np.float64)



def parse_args():
    parser = argparse.ArgumentParser(description='Tune the reference radius for the Monte Carlo score estimator.')
    parser.add_argument('--particles', type=int, default=10_000_000, help='Number of Monte Carlo particles in each realization.')
    parser.add_argument('--particle-sets', type=int, default=50, help='Number of independent particle realizations.')
    parser.add_argument('--seed', type=int, default=42, help='Random-number seed.')
    return parser.parse_args()


def load_centers():
    data_file = Path(__file__).resolve().parent / 'data' / 'initialData.mat'
    return np.asarray(scipy.io.loadmat(data_file)['initialData'], dtype=np.float64)


def standard_normal_density(points):
    return np.exp(-0.5 * np.sum(points**2, axis=1)) / (2.0 * np.pi)


def tune_reference_radius(centers, reference_radii, num_particles, num_particle_sets, seed):
    rng = np.random.default_rng(seed)
    exact_scores = -centers
    densities = standard_normal_density(centers)
    reference_density = np.median(densities)
    local_radii = reference_radii[:, None] * (reference_density / densities[None, :]) ** (1.0 / 8.0)

    squared_error_sums = np.zeros(reference_radii.size, dtype=np.float64)
    empty_ball_counts = np.zeros(reference_radii.size, dtype=np.int64)

    for particle_set in range(num_particle_sets):
        positions = rng.standard_normal((num_particles, 2), dtype=np.float64)

        for center_index, center in enumerate(centers):
            displacements = positions - center
            distances_squared = np.einsum('ij,ij->i', displacements, displacements)

            for radius_index in range(reference_radii.size):
                radius = local_radii[radius_index, center_index]
                inside = distances_squared < radius**2

                if np.any(inside):
                    estimated_score = 4.0 * np.mean(displacements[inside], axis=0) / radius**2
                else:
                    estimated_score = np.zeros(2, dtype=np.float64)
                    empty_ball_counts[radius_index] += 1

                score_error = estimated_score - exact_scores[center_index]
                squared_error_sums[radius_index] += np.dot(score_error, score_error)

        print(f'Completed particle set {particle_set + 1}/{num_particle_sets}')

    num_evaluations = num_particle_sets * centers.shape[0]
    mean_squared_errors = squared_error_sums / num_evaluations
    empty_ball_fractions = empty_ball_counts.astype(np.float64) / num_evaluations
    return mean_squared_errors, empty_ball_fractions, local_radii, reference_density


def main():
    args = parse_args()
    centers = load_centers()
    errors, empty_fractions, local_radii, reference_density = tune_reference_radius(
        centers, DEFAULT_RADII, args.particles, args.particle_sets, args.seed
    )
    optimal_index = int(np.argmin(errors))
    optimal_radius = DEFAULT_RADII[optimal_index]
    calibrated_constant = optimal_radius * (args.particles * reference_density) ** (1.0 / 8.0)

    print('\nReference radius    MSE              Empty-ball fraction    Local-radius range')
    for index, radius in enumerate(DEFAULT_RADII):
        radius_range = f'[{local_radii[index].min():.6g}, {local_radii[index].max():.6g}]'
        print(f'{radius:16.6g}    {errors[index]:.8e}    {empty_fractions[index]:19.6f}    {radius_range}')

    print(f'\nReference density: {reference_density:.8e}')
    print(f'Optimal reference radius: {optimal_radius:.8e}')
    print(f'Calibrated constant C: {calibrated_constant:.8e}')


if __name__ == '__main__':
    main()
