import argparse
from pathlib import Path

import numpy as np
import scipy.io


# STANDARD_REFERENCE_RADII = np.logspace(np.log10(0.25), np.log10(0.26), 10, dtype=np.float64)
STANDARD_REFERENCE_RADII = np.logspace(-2, 0, 10, dtype=np.float64)
SHIFTED_MEAN = np.ones(2, dtype=np.float64)
SHIFTED_VARIANCE = 2.0



def parse_args():
    parser = argparse.ArgumentParser(description='Tune the reference radius for the Monte Carlo score estimator.')
    parser.add_argument('--particles', type=int, default=10_000_000, help='Number of Monte Carlo particles in each realization.')
    parser.add_argument('--particle-sets', type=int, default=20, help='Number of independent particle realizations.')
    parser.add_argument('--seed', type=int, default=42, help='Random-number seed.')
    parser.add_argument('--setup', choices=('standard', 'shifted'), default='standard', help='Gaussian setup: standard normal or variance 2 with mean (1, 1).')
    return parser.parse_args()


def load_centers():
    data_file = Path(__file__).resolve().parent / 'data' / 'initialData.mat'
    return np.asarray(scipy.io.loadmat(data_file)['initialData'], dtype=np.float64)


def isotropic_gaussian_density(points, mean, variance):
    centered_points = points - mean
    return np.exp(-0.5 * np.sum(centered_points**2, axis=1) / variance) / (2.0 * np.pi * variance)


def tune_reference_radius(centers, reference_radii, num_particles, num_particle_sets, seed, mean, variance):
    rng = np.random.default_rng(seed)
    standard_deviation = np.sqrt(variance)
    exact_scores = -(centers - mean) / variance
    densities = isotropic_gaussian_density(centers, mean, variance)
    reference_density = np.median(densities)
    local_radii = reference_radii[:, None] * (reference_density / densities[None, :]) ** (1.0 / 8.0)

    squared_error_sums = np.zeros(reference_radii.size, dtype=np.float64)
    empty_ball_counts = np.zeros(reference_radii.size, dtype=np.int64)

    for particle_set in range(num_particle_sets):
        positions = mean + standard_deviation * rng.standard_normal((num_particles, 2), dtype=np.float64)

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
    if args.setup == 'standard':
        mean = np.zeros(2, dtype=np.float64)
        variance = 1.0
    else:
        mean = SHIFTED_MEAN
        variance = SHIFTED_VARIANCE

    standard_deviation = np.sqrt(variance)
    centers = mean + standard_deviation * load_centers()
    reference_radii = standard_deviation * STANDARD_REFERENCE_RADII
    errors, empty_fractions, local_radii, reference_density = tune_reference_radius(
        centers, reference_radii, args.particles, args.particle_sets, args.seed, mean, variance
    )
    optimal_index = int(np.argmin(errors))
    optimal_radius = reference_radii[optimal_index]
    constant_scale = variance ** (3.0 / 8.0)
    coordinate_constant = optimal_radius * (args.particles * reference_density) ** (1.0 / 8.0)
    calibrated_constant = coordinate_constant / constant_scale

    print(f'Gaussian setup: {args.setup}')
    print(f'Gaussian mean: {mean}')
    print(f'Gaussian variance: {variance:.8e}')
    print(f'Gaussian standard deviation: {standard_deviation:.8e}')
    print('\nReference radius    MSE              Empty-ball fraction    Local-radius range')
    for index, radius in enumerate(reference_radii):
        radius_range = f'[{local_radii[index].min():.6g}, {local_radii[index].max():.6g}]'
        print(f'{radius:16.6g}    {errors[index]:.8e}    {empty_fractions[index]:19.6f}    {radius_range}')

    print(f'\nReference density: {reference_density:.8e}')
    print(f'Optimal reference radius: {optimal_radius:.8e}')
    print(f'Coordinate-dependent constant: {coordinate_constant:.8e}')
    print(f'Calibrated scale-normalized constant C: {calibrated_constant:.8e}')

    if args.setup == 'shifted':
        density_scale = variance**-1.0
        radius_scale = standard_deviation

        print('\nTheoretical scaling relative to the standard setup')
        print(f'rho_ref scale (shifted/standard): {density_scale:.8e}')
        print(f'r_ref scale (shifted/standard): {radius_scale:.8e}')
        print(f'C scale (shifted/standard): {constant_scale:.8e}')
        print('\nShifted results mapped back to the standard scale')
        print(f'Scaled reference density: {reference_density / density_scale:.8e}')
        print(f'Scaled optimal reference radius: {optimal_radius / radius_scale:.8e}')
        print(f'Scaled coordinate-dependent constant: {coordinate_constant / constant_scale:.8e}')


if __name__ == '__main__':
    main()
