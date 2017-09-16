#ifndef RESAMPLERS
#define RESAMPLERS

#include "PF_data.h"
#include "particles.h"
#include "../sample_funcs.h"
#include "PF_utils.h"
#include "../utils.h"


#define MAX(a,b) (((a)>(b))?(a):(b))

/* base class for common functions*/

template <arma::uvec (*sample_func)(const arma::uword, arma::vec&)>
class resampler_base {
protected:
  inline static arma::uvec sample(
      const PF_data &data, arma::vec &probs, const double ESS, bool &did_resample){
    if(probs.n_elem != data.N_fw_n_bw){
      if(data.debug > 1){
        data.log(2) << "Subsampling " << probs.n_elem << " to get "
                    << data.N_fw_n_bw <<  " using re-sampling weights. ESS of re-sampling weights are "
                    << ESS;
      }

      did_resample = true;
      return(sample_func(data.N_fw_n_bw, probs));
    }

    if(ESS < data.forward_backward_ESS_threshold){

      if(data.debug > 1){
        data.log(2) << "ESS of re-sampling weights is below threshold (" << ESS << " < "
                    << data.forward_backward_ESS_threshold << "). Re-sampling";
      }

      if(data.debug > 2){
        data.log(3) << "Re-sampling " << data.N_fw_n_bw << " indices "
                    << " from " <<  probs.n_elem << " elements "
                    << " with " <<  arma::max(probs) << " as the higest probability";
      }

      did_resample = true;
      return(sample_func(data.N_fw_n_bw, probs));
    }

    if(data.debug > 1){
      data.log(2) << "ESS of re-sampling weights is greater than threshold (" << ESS << " >= "
                  << data.forward_backward_ESS_threshold << "). No re-sampling needed";
    }

    did_resample = false;
    return(arma::linspace<arma::uvec>(0, data.N_fw_n_bw - 1, data.N_fw_n_bw));
  }
};

/*
 Each "resampler" has a static function which:
    1) computes log re-sampling weights assuming that weights are computed
    2) samples according to re-sampling weights
    3) returns re-sampled indices if effective sample size is low. Otherwise
       return an index for each element in the cloud. A boolean argument is
       used to indicate if sampling is made
*/

/*
 Non-auxiliary particle filter where we just to the weights as they are. I.e.
 set:
  beta_j = w_{j - 1}
*/

template<typename densities, bool is_forward>
class None_AUX_resampler : private resampler_base<systematic_resampling> {
public:
  inline static nothing resampler(
      const PF_data &data, cloud &PF_cloud, unsigned int t, arma::uvec &outcome,
      bool &did_resample){
    /* Compute effective sample size (ESS) */
    arma::vec weights(PF_cloud.size());
    double ESS = 0;
    auto w = weights.begin();
    for(auto it = PF_cloud.begin(); it != PF_cloud.end(); ++it, ++w){
      /* No need to update weights */
      it->log_resampling_weight = it->log_weight;

      *w = exp(it->log_resampling_weight);
      ESS += *w * *w;
    }
    ESS = 1/ ESS;

    outcome = sample(data, weights, ESS, did_resample);

    return nothing();
  }
};

/*
  Auxiliary particle filter with weights as in the end of page 462 of:
    Fearnhead, P., Wyncoll, D., & Tawn, J. (2010). A sequential smoothing algorithm with linear computational cost. Biometrika, 97(2), 447-464.
*/

template<typename densities, bool is_forward>
class AUX_resampler_normal_approx_w_cloud_mean : private resampler_base<systematic_resampling> {
public:
  inline static input_for_normal_apprx_w_cloud_mean resampler(
      const PF_data &data, cloud &PF_cloud, unsigned int t, arma::uvec &outcome,
      bool &did_resample){
    /* Find weighted mean estimate */
    arma::vec alpha_bar = PF_cloud.get_weigthed_mean();

    /* compute means and covariances */
    auto &Q = data.Q_proposal;
    auto ans = compute_mu_n_Sigma_from_normal_apprx_w_cloud_mean
      <densities, is_forward>
      (data, t, Q, alpha_bar, PF_cloud);

    const arma::mat *Q_numerator_chol_inv;
    arma::vec *a_0_scaled;
    double bw_w1, bw_w2;
    if(is_forward){
      Q_numerator_chol_inv = &data.Q.chol_inv;

    } else {
      bw_w1 = (double)(t + 1) / (t + 2);
      bw_w2 =              1. / (t + 2);
      arma::vec tmp_vec = data.a_0 * bw_w2;

      a_0_scaled = new arma::vec(std::move(tmp_vec));
      arma::mat tmp = data.Q.chol_inv / sqrt((double)(t + 1) / (t + 2));
      Q_numerator_chol_inv = new arma::mat(std::move(tmp));

    }

    /* Compute sampling weights*/
    double max_weight =  -std::numeric_limits<double>::max();
    auto it_cl = PF_cloud.begin();
    auto it_mu_j = ans.mu_js.begin();
    unsigned int n_elem = PF_cloud.size();
    for(unsigned int i = 0; i != n_elem; ++i, ++it_cl, ++it_mu_j){
      double log_prob_y_given_state = densities::log_prob_y_given_state(
        data, *it_mu_j, t);
      double log_prop_transition;
      if(is_forward){
        log_prop_transition = dmvnrm_log(
          *it_mu_j, it_cl->get_state(), *Q_numerator_chol_inv);

      } else {
        arma::vec mean = bw_w1 * it_cl->get_state() + *a_0_scaled;
        log_prop_transition = dmvnrm_log(
          *it_mu_j, mean, *Q_numerator_chol_inv);

      }

      double log_prop_proposal = dmvnrm_log(
        *it_mu_j, *it_mu_j /* Notice same */, ans.sigma_chol_inv /* Notice Sigma*/);

      it_cl->log_resampling_weight =
        it_cl->log_weight + log_prop_transition + log_prob_y_given_state
        - log_prop_proposal;

      max_weight = MAX(it_cl->log_resampling_weight, max_weight);
    }

    if(!is_forward){
      delete Q_numerator_chol_inv;
      delete a_0_scaled;

    }

    auto norm_out = normalize_log_resampling_weight<true, true>(PF_cloud, max_weight);
    outcome = sample(data, norm_out.weights, norm_out.ESS, did_resample);

    return ans;
  }
};

/*
  Auxiliary particle filter with weights as in the end of page 462 of the
  following paper with Taylor expansion around the parent particle:
    Fearnhead, P., Wyncoll, D., & Tawn, J. (2010). A sequential smoothing algorithm with linear computational cost. Biometrika, 97(2), 447-464.
*/

template<typename densities, bool is_forward>
class AUX_resampler_normal_approx_w_particles : private resampler_base<systematic_resampling> {
public:
  inline static input_for_normal_apprx_w_particle_mean resampler(
      const PF_data &data, cloud &PF_cloud, unsigned int t, arma::uvec &outcome,
      bool &did_resample){
    /* compute means and covariances */
    auto &Q = data.Q_proposal;
    auto ans = compute_mu_n_Sigma_from_normal_apprx_w_particles
      <densities, is_forward>
      (data, t, Q, PF_cloud);

    arma::vec *a_0_scaled;
    double bw_w1, bw_w2;
    const arma::mat *Q_numerator_chol_inv;
    if(is_forward){
      Q_numerator_chol_inv = &data.Q.chol_inv;

    } else {
      bw_w1 = (double)(t + 1) / (t + 2);
      bw_w2 =              1. / (t + 2);
      arma::vec tmp_vec = data.a_0 * bw_w2;

      a_0_scaled = new arma::vec(std::move(tmp_vec));
      arma::mat tmp = data.Q.chol_inv / sqrt((double)(t + 1) / (t + 2));
      Q_numerator_chol_inv = new arma::mat(std::move(tmp));

    }

    /* Compute sampling weights */
    double max_weight =  -std::numeric_limits<double>::max();

    unsigned int n_elem = PF_cloud.size();
    arma::uvec r_set = get_risk_set(data, t);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for(unsigned int i = 0; i < n_elem; ++i){
      auto it_cl = &PF_cloud[i];
      auto it_ans = &ans[i];

      double log_prob_y_given_state = densities::log_prob_y_given_state(
        data, it_ans->mu, t, r_set, false);

      double log_prop_transition;
      if(is_forward){
        log_prop_transition = dmvnrm_log(
          it_ans->mu, it_cl->get_state(), *Q_numerator_chol_inv);

      } else {
        arma::vec mean = bw_w1 * it_cl->get_state() + *a_0_scaled;
        log_prop_transition = dmvnrm_log(
          it_ans->mu, mean, *Q_numerator_chol_inv);

      }

      double log_prop_proposal = dmvnrm_log(
        it_ans->mu, it_ans->mu /* Notice same */, it_ans->sigma_chol_inv /* Notice Sigma*/);

      it_cl->log_resampling_weight =
      it_cl->log_weight + log_prop_transition + log_prob_y_given_state
        - log_prop_proposal;
#ifdef _OPENMP
#pragma omp critical(resampler_aux_norm_particle_lock)
{
#endif
      max_weight = MAX(it_cl->log_resampling_weight, max_weight);
#ifdef _OPENMP
}
#endif
    }

    if(!is_forward){
      delete Q_numerator_chol_inv;
      delete a_0_scaled;

    }

    auto norm_out = normalize_log_resampling_weight<true, true>(PF_cloud, max_weight);
    outcome = sample(data, norm_out.weights, norm_out.ESS, did_resample);

    return ans;
  }
};



#undef MAX
#endif