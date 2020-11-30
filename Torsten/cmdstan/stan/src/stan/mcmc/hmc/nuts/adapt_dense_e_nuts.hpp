#ifndef STAN_MCMC_HMC_NUTS_ADAPT_DENSE_E_NUTS_HPP
#define STAN_MCMC_HMC_NUTS_ADAPT_DENSE_E_NUTS_HPP

#include <stan/callbacks/logger.hpp>
#include <stan/mcmc/stepsize_covar_adapter.hpp>
#include <stan/mcmc/hmc/nuts/dense_e_nuts.hpp>
#include <stan/mcmc/cross_chain/mpi_cross_chain_adapter.hpp>

namespace stan {
namespace mcmc {
/**
 * The No-U-Turn sampler (NUTS) with multinomial sampling
 * with a Gaussian-Euclidean disintegration and adaptive
 * dense metric and adaptive step size
 */
template <class Model, class BaseRNG>
class adapt_dense_e_nuts : public dense_e_nuts<Model, BaseRNG>,
                           public stepsize_covar_adapter,
                           public mpi_cross_chain_adapter {
 public:
  adapt_dense_e_nuts(const Model& model, BaseRNG& rng)
      : dense_e_nuts<Model, BaseRNG>(model, rng),
        stepsize_covar_adapter(model.num_params_r()) {}

  ~adapt_dense_e_nuts() {}

  sample transition(sample& init_sample, callbacks::logger& logger) {
    sample s = dense_e_nuts<Model, BaseRNG>::transition(init_sample, logger);

    if (this->adapt_flag_) {
      this->stepsize_adaptation_.learn_stepsize(this->nom_epsilon_,
                                                s.accept_stat());

      if (this -> use_cross_chain_adapt()) {
        /// cross chain adapter has its own var adaptor so needs to add sample
        this -> add_cross_chain_sample(s.log_prob(), this -> z().q);
        bool update = this -> cross_chain_adaptation(this -> z().inv_e_metric_, logger);
        if (update) {
          // this->init_stepsize(logger);
          double new_stepsize = this -> cross_chain_stepsize(this -> get_nominal_stepsize());
          this -> set_nominal_stepsize(new_stepsize);
          this->stepsize_adaptation_.set_mu(log(10 * this->nom_epsilon_));
          this->stepsize_adaptation_.restart();
        }
      } else {
        bool update = this->covar_adaptation_.learn_covariance(this->z_.inv_e_metric_,
                                                               this->z_.q);

        if (update) {
          this->init_stepsize(logger);
          this->stepsize_adaptation_.set_mu(log(10 * this->nom_epsilon_));
          this->stepsize_adaptation_.restart();
        }
      }
    }
    return s;
  }

  void disengage_adaptation() {
    base_adapter::disengage_adaptation();
    this->stepsize_adaptation_.complete_adaptation(this->nom_epsilon_);
  }
};

}  // namespace mcmc
}  // namespace stan
#endif