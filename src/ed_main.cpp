#define main ftlm_lanczos_driver_main
#include "main.cpp"
#undef main

namespace {

void accumulate_exact_block(
    const std::vector<double>& eigenvalues,
    int particles,
    double beta,
    const std::vector<double>& mu_values,
    std::vector<double>& emin_by_mu,
    std::vector<double>& z_scaled,
    std::vector<double>& n_scaled,
    std::vector<double>& n2_scaled) {
  for (std::size_t imu = 0; imu < mu_values.size(); ++imu) {
    const double mu = mu_values[imu];
    double block_min = std::numeric_limits<double>::infinity();
    for (double eigenvalue : eigenvalues) {
      block_min = std::min(block_min, eigenvalue - mu * static_cast<double>(particles));
    }
    if (!std::isfinite(block_min)) {
      continue;
    }

    if (!std::isfinite(emin_by_mu[imu])) {
      emin_by_mu[imu] = block_min;
    } else if (block_min < emin_by_mu[imu]) {
      const double scale = std::exp(-beta * (emin_by_mu[imu] - block_min));
      z_scaled[imu] *= scale;
      n_scaled[imu] *= scale;
      n2_scaled[imu] *= scale;
      emin_by_mu[imu] = block_min;
    }

    for (double eigenvalue : eigenvalues) {
      const double shifted = eigenvalue - mu * static_cast<double>(particles) - emin_by_mu[imu];
      const double weight = std::exp(-beta * shifted);
      const double particle_count = static_cast<double>(particles);
      z_scaled[imu] += weight;
      n_scaled[imu] += particle_count * weight;
      n2_scaled[imu] += particle_count * particle_count * weight;
    }
  }
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const Params params = parse_args(argc, argv);
    const Lattice lattice = build_lattice(params.lx, params.ly);
    const std::vector<double> mu_values =
        linspace(params.mu_min, params.mu_max, params.mu_count);

    std::cout << "Exact diagonalization using the FTLM sector and k-block builder\n";
    std::cout << "Running ED on a " << params.lx << "x" << params.ly
              << " rectangular lattice (" << lattice.sites << " sites)\n";
    std::cout << "beta=" << params.beta << " U=" << params.u << " V=" << params.v
              << " tx=" << params.tx << " ty=" << params.ty << " tp=" << params.tp
              << " phix=" << params.phix << " phiy=" << params.phiy << "\n";

    std::vector<double> densities(mu_values.size(), 0.0);
    std::vector<double> charge_correlations(mu_values.size(), 0.0);
    std::vector<double> compressibilities(mu_values.size(), 0.0);
    std::vector<double> partitions(mu_values.size(), 0.0);
    std::vector<double> log_partitions(mu_values.size(), 0.0);
    std::vector<double> emin_by_mu(
        mu_values.size(), std::numeric_limits<double>::infinity());
    std::vector<double> z_scaled(mu_values.size(), 0.0);
    std::vector<double> n_scaled(mu_values.size(), 0.0);
    std::vector<double> n2_scaled(mu_values.size(), 0.0);

    for (int n_up = 0; n_up <= lattice.sites; ++n_up) {
      for (int n_down = 0; n_down <= lattice.sites; ++n_down) {
        const SectorBasis sector = build_sector_basis(lattice, n_up, n_down, params);
        if (sector.full_dim == 0) {
          continue;
        }
        if (sector.full_dim > params.max_sector_dim) {
          std::ostringstream msg;
          msg << "Sector (Nup=" << n_up << ", Ndown=" << n_down << ") has full dimension "
              << sector.full_dim << ", which exceeds --max-sector-dim=" << params.max_sector_dim;
          die(msg.str());
        }

        std::size_t active_sum = 0;
        for (int my = 0; my < lattice.ly; ++my) {
          for (int mx = 0; mx < lattice.lx; ++mx) {
            const MomentumBlock block = build_momentum_block(sector, lattice, mx, my);
            active_sum += block.basis_dim;
            if (block.basis_dim == 0) {
              continue;
            }
            const std::vector<double> eigenvalues = diagonalize_block_exact(block);
            accumulate_exact_block(
                eigenvalues,
                block.particles,
                params.beta,
                mu_values,
                emin_by_mu,
                z_scaled,
                n_scaled,
                n2_scaled);
          }
        }

        std::cout << "  sector Nup=" << n_up << " Ndown=" << n_down
                  << " full_dim=" << sector.full_dim
                  << " parents=" << sector.parents.size()
                  << " active_k_dim_sum=" << active_sum << "\n";
      }
    }

    for (std::size_t imu = 0; imu < mu_values.size(); ++imu) {
      partitions[imu] = z_scaled[imu];
      log_partitions[imu] =
          std::log(z_scaled[imu]) - params.beta * emin_by_mu[imu];
      const double sites = static_cast<double>(lattice.sites);
      const double mean_particles = n_scaled[imu] / z_scaled[imu];
      const double mean_particles_squared = n2_scaled[imu] / z_scaled[imu];
      const double charge_correlation =
          (mean_particles_squared - mean_particles * mean_particles) / sites;
      densities[imu] = mean_particles / sites;
      charge_correlations[imu] = std::max(0.0, charge_correlation);
      compressibilities[imu] = params.beta * charge_correlations[imu];
      std::cout << "mu=" << std::setw(12) << mu_values[imu] << "  n=" << std::setw(12)
                << densities[imu] << "  kappa=" << std::setw(12)
                << compressibilities[imu] << "\n";
    }

    write_results(
        params.output,
        mu_values,
        densities,
        charge_correlations,
        compressibilities,
        partitions,
        log_partitions);
    std::cout << "Wrote " << params.output << "\n";
    return 0;
  } catch (const std::exception& ex) {
    std::cerr << "Error: " << ex.what() << "\n";
    return 1;
  }
}
