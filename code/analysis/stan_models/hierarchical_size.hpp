
// Code generated by stanc v2.30.1
#include <stan/model/model_header.hpp>
namespace hierarchical_size_model_namespace {

using stan::model::model_base_crtp;
using namespace stan::math;


stan::math::profile_map profiles__;
static constexpr std::array<const char*, 35> locations_array__ = 
{" (found before start of program)",
 " (in '/Users/gchure/Dropbox/git/postdoc_projects/size_control/code/analysis/stan_models/hierarchical_size.stan', line 19, column 4 to column 20)",
 " (in '/Users/gchure/Dropbox/git/postdoc_projects/size_control/code/analysis/stan_models/hierarchical_size.stan', line 20, column 4 to column 21)",
 " (in '/Users/gchure/Dropbox/git/postdoc_projects/size_control/code/analysis/stan_models/hierarchical_size.stan', line 21, column 4 to column 32)",
 " (in '/Users/gchure/Dropbox/git/postdoc_projects/size_control/code/analysis/stan_models/hierarchical_size.stan', line 22, column 4 to column 33)",
 " (in '/Users/gchure/Dropbox/git/postdoc_projects/size_control/code/analysis/stan_models/hierarchical_size.stan', line 25, column 4 to column 23)",
 " (in '/Users/gchure/Dropbox/git/postdoc_projects/size_control/code/analysis/stan_models/hierarchical_size.stan', line 26, column 4 to column 24)",
 " (in '/Users/gchure/Dropbox/git/postdoc_projects/size_control/code/analysis/stan_models/hierarchical_size.stan', line 29, column 4 to column 35)",
 " (in '/Users/gchure/Dropbox/git/postdoc_projects/size_control/code/analysis/stan_models/hierarchical_size.stan', line 30, column 4 to column 36)",
 " (in '/Users/gchure/Dropbox/git/postdoc_projects/size_control/code/analysis/stan_models/hierarchical_size.stan', line 35, column 4 to column 30)",
 " (in '/Users/gchure/Dropbox/git/postdoc_projects/size_control/code/analysis/stan_models/hierarchical_size.stan', line 36, column 4 to column 33)",
 " (in '/Users/gchure/Dropbox/git/postdoc_projects/size_control/code/analysis/stan_models/hierarchical_size.stan', line 37, column 4 to column 31)",
 " (in '/Users/gchure/Dropbox/git/postdoc_projects/size_control/code/analysis/stan_models/hierarchical_size.stan', line 38, column 4 to column 34)",
 " (in '/Users/gchure/Dropbox/git/postdoc_projects/size_control/code/analysis/stan_models/hierarchical_size.stan', line 41, column 4 to column 49)",
 " (in '/Users/gchure/Dropbox/git/postdoc_projects/size_control/code/analysis/stan_models/hierarchical_size.stan', line 42, column 4 to column 52)",
 " (in '/Users/gchure/Dropbox/git/postdoc_projects/size_control/code/analysis/stan_models/hierarchical_size.stan', line 43, column 4 to column 31)",
 " (in '/Users/gchure/Dropbox/git/postdoc_projects/size_control/code/analysis/stan_models/hierarchical_size.stan', line 44, column 4 to column 32)",
 " (in '/Users/gchure/Dropbox/git/postdoc_projects/size_control/code/analysis/stan_models/hierarchical_size.stan', line 46, column 4 to column 64)",
 " (in '/Users/gchure/Dropbox/git/postdoc_projects/size_control/code/analysis/stan_models/hierarchical_size.stan', line 47, column 4 to column 67)",
 " (in '/Users/gchure/Dropbox/git/postdoc_projects/size_control/code/analysis/stan_models/hierarchical_size.stan', line 3, column 4 to column 19)",
 " (in '/Users/gchure/Dropbox/git/postdoc_projects/size_control/code/analysis/stan_models/hierarchical_size.stan', line 4, column 4 to column 19)",
 " (in '/Users/gchure/Dropbox/git/postdoc_projects/size_control/code/analysis/stan_models/hierarchical_size.stan', line 5, column 30 to column 31)",
 " (in '/Users/gchure/Dropbox/git/postdoc_projects/size_control/code/analysis/stan_models/hierarchical_size.stan', line 5, column 4 to column 33)",
 " (in '/Users/gchure/Dropbox/git/postdoc_projects/size_control/code/analysis/stan_models/hierarchical_size.stan', line 8, column 20 to column 21)",
 " (in '/Users/gchure/Dropbox/git/postdoc_projects/size_control/code/analysis/stan_models/hierarchical_size.stan', line 8, column 4 to column 30)",
 " (in '/Users/gchure/Dropbox/git/postdoc_projects/size_control/code/analysis/stan_models/hierarchical_size.stan', line 9, column 20 to column 21)",
 " (in '/Users/gchure/Dropbox/git/postdoc_projects/size_control/code/analysis/stan_models/hierarchical_size.stan', line 9, column 4 to column 31)",
 " (in '/Users/gchure/Dropbox/git/postdoc_projects/size_control/code/analysis/stan_models/hierarchical_size.stan', line 13, column 11 to column 12)",
 " (in '/Users/gchure/Dropbox/git/postdoc_projects/size_control/code/analysis/stan_models/hierarchical_size.stan', line 13, column 4 to column 56)",
 " (in '/Users/gchure/Dropbox/git/postdoc_projects/size_control/code/analysis/stan_models/hierarchical_size.stan', line 14, column 11 to column 12)",
 " (in '/Users/gchure/Dropbox/git/postdoc_projects/size_control/code/analysis/stan_models/hierarchical_size.stan', line 14, column 4 to column 59)",
 " (in '/Users/gchure/Dropbox/git/postdoc_projects/size_control/code/analysis/stan_models/hierarchical_size.stan', line 25, column 11 to column 12)",
 " (in '/Users/gchure/Dropbox/git/postdoc_projects/size_control/code/analysis/stan_models/hierarchical_size.stan', line 26, column 11 to column 12)",
 " (in '/Users/gchure/Dropbox/git/postdoc_projects/size_control/code/analysis/stan_models/hierarchical_size.stan', line 29, column 20 to column 21)",
 " (in '/Users/gchure/Dropbox/git/postdoc_projects/size_control/code/analysis/stan_models/hierarchical_size.stan', line 30, column 20 to column 21)"};




class hierarchical_size_model final : public model_base_crtp<hierarchical_size_model> {

 private:
  int J;
  int N;
  std::vector<int> idx;
  Eigen::Matrix<double, -1, 1> widths_data__;
  Eigen::Matrix<double, -1, 1> lengths_data__;
  Eigen::Matrix<double, -1, 1> widths_uncentered_data__;
  Eigen::Matrix<double, -1, 1> lengths_uncentered_data__; 
  Eigen::Map<Eigen::Matrix<double, -1, 1>> widths{nullptr, 0};
  Eigen::Map<Eigen::Matrix<double, -1, 1>> lengths{nullptr, 0};
  Eigen::Map<Eigen::Matrix<double, -1, 1>> widths_uncentered{nullptr, 0};
  Eigen::Map<Eigen::Matrix<double, -1, 1>> lengths_uncentered{nullptr, 0};
 
 public:
  ~hierarchical_size_model() { }
  
  inline std::string model_name() const final { return "hierarchical_size_model"; }

  inline std::vector<std::string> model_compile_info() const noexcept {
    return std::vector<std::string>{"stanc_version = stanc3 v2.30.1", "stancflags = "};
  }
  
  
  hierarchical_size_model(stan::io::var_context& context__,
                          unsigned int random_seed__ = 0,
                          std::ostream* pstream__ = nullptr) : model_base_crtp(0) {
    int current_statement__ = 0;
    using local_scalar_t__ = double ;
    boost::ecuyer1988 base_rng__ = 
        stan::services::util::create_rng(random_seed__, 0);
    (void) base_rng__;  // suppress unused var warning
    static constexpr const char* function__ = "hierarchical_size_model_namespace::hierarchical_size_model";
    (void) function__;  // suppress unused var warning
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    (void) DUMMY_VAR__;  // suppress unused var warning
    try {
      int pos__ = std::numeric_limits<int>::min();
      pos__ = 1;
      current_statement__ = 19;
      context__.validate_dims("data initialization","J","int",
           std::vector<size_t>{});
      J = std::numeric_limits<int>::min();
      
      
      current_statement__ = 19;
      J = context__.vals_i("J")[(1 - 1)];
      current_statement__ = 19;
      stan::math::check_greater_or_equal(function__, "J", J, 1);
      current_statement__ = 20;
      context__.validate_dims("data initialization","N","int",
           std::vector<size_t>{});
      N = std::numeric_limits<int>::min();
      
      
      current_statement__ = 20;
      N = context__.vals_i("N")[(1 - 1)];
      current_statement__ = 20;
      stan::math::check_greater_or_equal(function__, "N", N, 1);
      current_statement__ = 21;
      stan::math::validate_non_negative_index("idx", "N", N);
      current_statement__ = 22;
      context__.validate_dims("data initialization","idx","int",
           std::vector<size_t>{static_cast<size_t>(N)});
      idx = std::vector<int>(N, std::numeric_limits<int>::min());
      
      
      current_statement__ = 22;
      idx = context__.vals_i("idx");
      current_statement__ = 22;
      stan::math::check_greater_or_equal(function__, "idx", idx, 1);
      current_statement__ = 22;
      stan::math::check_less_or_equal(function__, "idx", idx, J);
      current_statement__ = 23;
      stan::math::validate_non_negative_index("widths", "N", N);
      current_statement__ = 24;
      context__.validate_dims("data initialization","widths","double",
           std::vector<size_t>{static_cast<size_t>(N)});
      widths_data__ = 
        Eigen::Matrix<double, -1, 1>::Constant(N,
          std::numeric_limits<double>::quiet_NaN());
      new (&widths) Eigen::Map<Eigen::Matrix<double, -1, 1>>(widths_data__.data(), N);
        
      
      {
        std::vector<local_scalar_t__> widths_flat__;
        current_statement__ = 24;
        widths_flat__ = context__.vals_r("widths");
        current_statement__ = 24;
        pos__ = 1;
        current_statement__ = 24;
        for (int sym1__ = 1; sym1__ <= N; ++sym1__) {
          current_statement__ = 24;
          stan::model::assign(widths, widths_flat__[(pos__ - 1)],
            "assigning variable widths", stan::model::index_uni(sym1__));
          current_statement__ = 24;
          pos__ = (pos__ + 1);
        }
      }
      current_statement__ = 24;
      stan::math::check_greater_or_equal(function__, "widths", widths, 0);
      current_statement__ = 25;
      stan::math::validate_non_negative_index("lengths", "N", N);
      current_statement__ = 26;
      context__.validate_dims("data initialization","lengths","double",
           std::vector<size_t>{static_cast<size_t>(N)});
      lengths_data__ = 
        Eigen::Matrix<double, -1, 1>::Constant(N,
          std::numeric_limits<double>::quiet_NaN());
      new (&lengths) Eigen::Map<Eigen::Matrix<double, -1, 1>>(lengths_data__.data(), N);
        
      
      {
        std::vector<local_scalar_t__> lengths_flat__;
        current_statement__ = 26;
        lengths_flat__ = context__.vals_r("lengths");
        current_statement__ = 26;
        pos__ = 1;
        current_statement__ = 26;
        for (int sym1__ = 1; sym1__ <= N; ++sym1__) {
          current_statement__ = 26;
          stan::model::assign(lengths, lengths_flat__[(pos__ - 1)],
            "assigning variable lengths", stan::model::index_uni(sym1__));
          current_statement__ = 26;
          pos__ = (pos__ + 1);
        }
      }
      current_statement__ = 26;
      stan::math::check_greater_or_equal(function__, "lengths", lengths, 0);
      current_statement__ = 27;
      stan::math::validate_non_negative_index("widths_uncentered", "N", N);
      current_statement__ = 28;
      widths_uncentered_data__ = 
        Eigen::Matrix<double, -1, 1>::Constant(N,
          std::numeric_limits<double>::quiet_NaN());
      new (&widths_uncentered) Eigen::Map<Eigen::Matrix<double, -1, 1>>(widths_uncentered_data__.data(), N);
        
      
      current_statement__ = 28;
      stan::model::assign(widths_uncentered,
        stan::math::subtract(widths, stan::math::mean(widths)),
        "assigning variable widths_uncentered");
      current_statement__ = 29;
      stan::math::validate_non_negative_index("lengths_uncentered", "N", N);
      current_statement__ = 30;
      lengths_uncentered_data__ = 
        Eigen::Matrix<double, -1, 1>::Constant(N,
          std::numeric_limits<double>::quiet_NaN());
      new (&lengths_uncentered) Eigen::Map<Eigen::Matrix<double, -1, 1>>(lengths_uncentered_data__.data(), N);
        
      
      current_statement__ = 30;
      stan::model::assign(lengths_uncentered,
        stan::math::subtract(lengths, stan::math::mean(lengths)),
        "assigning variable lengths_uncentered");
      current_statement__ = 31;
      stan::math::validate_non_negative_index("width_mu", "J", J);
      current_statement__ = 32;
      stan::math::validate_non_negative_index("length_mu", "J", J);
      current_statement__ = 33;
      stan::math::validate_non_negative_index("width_sigma", "J", J);
      current_statement__ = 34;
      stan::math::validate_non_negative_index("length_sigma", "J", J);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
    num_params_r__ = 1 + 1 + 1 + 1 + J + J + J + J;
    
  }
  
  template <bool propto__, bool jacobian__ , typename VecR, typename VecI, 
  stan::require_vector_like_t<VecR>* = nullptr, 
  stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr> 
  inline stan::scalar_type_t<VecR> log_prob_impl(VecR& params_r__,
                                                 VecI& params_i__,
                                                 std::ostream* pstream__ = nullptr) const {
    using T__ = stan::scalar_type_t<VecR>;
    using local_scalar_t__ = T__;
    T__ lp__(0.0);
    stan::math::accumulator<T__> lp_accum__;
    stan::io::deserializer<local_scalar_t__> in__(params_r__, params_i__);
    int current_statement__ = 0;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    (void) DUMMY_VAR__;  // suppress unused var warning
    static constexpr const char* function__ = "hierarchical_size_model_namespace::log_prob";
    (void) function__;  // suppress unused var warning
    
    try {
      local_scalar_t__ width_mu_0 = DUMMY_VAR__;
      current_statement__ = 1;
      width_mu_0 = in__.template read<local_scalar_t__>();
      local_scalar_t__ length_mu_0 = DUMMY_VAR__;
      current_statement__ = 2;
      length_mu_0 = in__.template read<local_scalar_t__>();
      local_scalar_t__ width_sigma_0 = DUMMY_VAR__;
      current_statement__ = 3;
      width_sigma_0 = in__.template read_constrain_lb<local_scalar_t__, 
                        jacobian__>(0, lp__);
      local_scalar_t__ length_sigma_0 = DUMMY_VAR__;
      current_statement__ = 4;
      length_sigma_0 = in__.template read_constrain_lb<local_scalar_t__, 
                         jacobian__>(0, lp__);
      Eigen::Matrix<local_scalar_t__, -1, 1> width_mu =
         Eigen::Matrix<local_scalar_t__, -1, 1>::Constant(J, DUMMY_VAR__);
      current_statement__ = 5;
      width_mu = in__.template read<Eigen::Matrix<local_scalar_t__, -1, 1>>(
                   J);
      Eigen::Matrix<local_scalar_t__, -1, 1> length_mu =
         Eigen::Matrix<local_scalar_t__, -1, 1>::Constant(J, DUMMY_VAR__);
      current_statement__ = 6;
      length_mu = in__.template read<Eigen::Matrix<local_scalar_t__, -1, 1>>(
                    J);
      Eigen::Matrix<local_scalar_t__, -1, 1> width_sigma =
         Eigen::Matrix<local_scalar_t__, -1, 1>::Constant(J, DUMMY_VAR__);
      current_statement__ = 7;
      width_sigma = in__.template read_constrain_lb<
                      Eigen::Matrix<local_scalar_t__, -1, 1>, jacobian__>(0,
                      lp__, J);
      Eigen::Matrix<local_scalar_t__, -1, 1> length_sigma =
         Eigen::Matrix<local_scalar_t__, -1, 1>::Constant(J, DUMMY_VAR__);
      current_statement__ = 8;
      length_sigma = in__.template read_constrain_lb<
                       Eigen::Matrix<local_scalar_t__, -1, 1>, jacobian__>(0,
                       lp__, J);
      {
        current_statement__ = 9;
        lp_accum__.add(stan::math::normal_lpdf<propto__>(width_mu_0, 0, 1));
        current_statement__ = 10;
        lp_accum__.add(stan::math::normal_lpdf<propto__>(width_sigma_0, 0, 1));
        current_statement__ = 11;
        lp_accum__.add(stan::math::normal_lpdf<propto__>(length_mu_0, 0, 1));
        current_statement__ = 12;
        lp_accum__.add(
          stan::math::normal_lpdf<propto__>(length_sigma_0, 0, 1));
        current_statement__ = 13;
        lp_accum__.add(
          stan::math::normal_lpdf<propto__>(width_mu, width_mu_0,
            width_sigma_0));
        current_statement__ = 14;
        lp_accum__.add(
          stan::math::normal_lpdf<propto__>(length_mu, length_mu_0,
            length_sigma_0));
        current_statement__ = 15;
        lp_accum__.add(stan::math::normal_lpdf<propto__>(width_sigma, 0, 1));
        current_statement__ = 16;
        lp_accum__.add(stan::math::normal_lpdf<propto__>(length_sigma, 0, 1));
        current_statement__ = 17;
        lp_accum__.add(
          stan::math::normal_lpdf<propto__>(widths_uncentered,
            stan::model::rvalue(width_mu, "width_mu",
              stan::model::index_multi(idx)),
            stan::model::rvalue(width_sigma, "width_sigma",
              stan::model::index_multi(idx))));
        current_statement__ = 18;
        lp_accum__.add(
          stan::math::normal_lpdf<propto__>(lengths_uncentered,
            stan::model::rvalue(length_mu, "length_mu",
              stan::model::index_multi(idx)),
            stan::model::rvalue(length_sigma, "length_sigma",
              stan::model::index_multi(idx))));
      }
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
    lp_accum__.add(lp__);
    return lp_accum__.sum();
    } // log_prob_impl() 
    
  template <typename RNG, typename VecR, typename VecI, typename VecVar, 
  stan::require_vector_like_vt<std::is_floating_point, VecR>* = nullptr, 
  stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr, 
  stan::require_vector_vt<std::is_floating_point, VecVar>* = nullptr> 
  inline void write_array_impl(RNG& base_rng__, VecR& params_r__,
                               VecI& params_i__, VecVar& vars__,
                               const bool emit_transformed_parameters__ = true,
                               const bool emit_generated_quantities__ = true,
                               std::ostream* pstream__ = nullptr) const {
    using local_scalar_t__ = double;
    stan::io::deserializer<local_scalar_t__> in__(params_r__, params_i__);
    stan::io::serializer<local_scalar_t__> out__(vars__);
    static constexpr bool propto__ = true;
    (void) propto__;
    double lp__ = 0.0;
    (void) lp__;  // dummy to suppress unused var warning
    int current_statement__ = 0; 
    stan::math::accumulator<double> lp_accum__;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    constexpr bool jacobian__ = false;
    (void) DUMMY_VAR__;  // suppress unused var warning
    static constexpr const char* function__ = "hierarchical_size_model_namespace::write_array";
    (void) function__;  // suppress unused var warning
    
    try {
      double width_mu_0 = std::numeric_limits<double>::quiet_NaN();
      current_statement__ = 1;
      width_mu_0 = in__.template read<local_scalar_t__>();
      double length_mu_0 = std::numeric_limits<double>::quiet_NaN();
      current_statement__ = 2;
      length_mu_0 = in__.template read<local_scalar_t__>();
      double width_sigma_0 = std::numeric_limits<double>::quiet_NaN();
      current_statement__ = 3;
      width_sigma_0 = in__.template read_constrain_lb<local_scalar_t__, 
                        jacobian__>(0, lp__);
      double length_sigma_0 = std::numeric_limits<double>::quiet_NaN();
      current_statement__ = 4;
      length_sigma_0 = in__.template read_constrain_lb<local_scalar_t__, 
                         jacobian__>(0, lp__);
      Eigen::Matrix<double, -1, 1> width_mu =
         Eigen::Matrix<double, -1, 1>::Constant(J,
           std::numeric_limits<double>::quiet_NaN());
      current_statement__ = 5;
      width_mu = in__.template read<Eigen::Matrix<local_scalar_t__, -1, 1>>(
                   J);
      Eigen::Matrix<double, -1, 1> length_mu =
         Eigen::Matrix<double, -1, 1>::Constant(J,
           std::numeric_limits<double>::quiet_NaN());
      current_statement__ = 6;
      length_mu = in__.template read<Eigen::Matrix<local_scalar_t__, -1, 1>>(
                    J);
      Eigen::Matrix<double, -1, 1> width_sigma =
         Eigen::Matrix<double, -1, 1>::Constant(J,
           std::numeric_limits<double>::quiet_NaN());
      current_statement__ = 7;
      width_sigma = in__.template read_constrain_lb<
                      Eigen::Matrix<local_scalar_t__, -1, 1>, jacobian__>(0,
                      lp__, J);
      Eigen::Matrix<double, -1, 1> length_sigma =
         Eigen::Matrix<double, -1, 1>::Constant(J,
           std::numeric_limits<double>::quiet_NaN());
      current_statement__ = 8;
      length_sigma = in__.template read_constrain_lb<
                       Eigen::Matrix<local_scalar_t__, -1, 1>, jacobian__>(0,
                       lp__, J);
      out__.write(width_mu_0);
      out__.write(length_mu_0);
      out__.write(width_sigma_0);
      out__.write(length_sigma_0);
      out__.write(width_mu);
      out__.write(length_mu);
      out__.write(width_sigma);
      out__.write(length_sigma);
      if (stan::math::logical_negation((stan::math::primitive_value(
            emit_transformed_parameters__) || stan::math::primitive_value(
            emit_generated_quantities__)))) {
        return ;
      } 
      if (stan::math::logical_negation(emit_generated_quantities__)) {
        return ;
      } 
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
    } // write_array_impl() 
    
  template <typename VecVar, typename VecI, 
  stan::require_vector_t<VecVar>* = nullptr, 
  stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr> 
  inline void transform_inits_impl(VecVar& params_r__, VecI& params_i__,
                                   VecVar& vars__,
                                   std::ostream* pstream__ = nullptr) const {
    using local_scalar_t__ = double;
    stan::io::deserializer<local_scalar_t__> in__(params_r__, params_i__);
    stan::io::serializer<local_scalar_t__> out__(vars__);
    int current_statement__ = 0;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    
    try {
      int pos__ = std::numeric_limits<int>::min();
      pos__ = 1;
      local_scalar_t__ width_mu_0 = DUMMY_VAR__;
      width_mu_0 = in__.read<local_scalar_t__>();
      out__.write(width_mu_0);
      local_scalar_t__ length_mu_0 = DUMMY_VAR__;
      length_mu_0 = in__.read<local_scalar_t__>();
      out__.write(length_mu_0);
      local_scalar_t__ width_sigma_0 = DUMMY_VAR__;
      width_sigma_0 = in__.read<local_scalar_t__>();
      out__.write_free_lb(0, width_sigma_0);
      local_scalar_t__ length_sigma_0 = DUMMY_VAR__;
      length_sigma_0 = in__.read<local_scalar_t__>();
      out__.write_free_lb(0, length_sigma_0);
      Eigen::Matrix<local_scalar_t__, -1, 1> width_mu =
         Eigen::Matrix<local_scalar_t__, -1, 1>::Constant(J, DUMMY_VAR__);
      for (int sym1__ = 1; sym1__ <= J; ++sym1__) {
        stan::model::assign(width_mu, in__.read<local_scalar_t__>(),
          "assigning variable width_mu", stan::model::index_uni(sym1__));
      }
      out__.write(width_mu);
      Eigen::Matrix<local_scalar_t__, -1, 1> length_mu =
         Eigen::Matrix<local_scalar_t__, -1, 1>::Constant(J, DUMMY_VAR__);
      for (int sym1__ = 1; sym1__ <= J; ++sym1__) {
        stan::model::assign(length_mu, in__.read<local_scalar_t__>(),
          "assigning variable length_mu", stan::model::index_uni(sym1__));
      }
      out__.write(length_mu);
      Eigen::Matrix<local_scalar_t__, -1, 1> width_sigma =
         Eigen::Matrix<local_scalar_t__, -1, 1>::Constant(J, DUMMY_VAR__);
      for (int sym1__ = 1; sym1__ <= J; ++sym1__) {
        stan::model::assign(width_sigma, in__.read<local_scalar_t__>(),
          "assigning variable width_sigma", stan::model::index_uni(sym1__));
      }
      out__.write_free_lb(0, width_sigma);
      Eigen::Matrix<local_scalar_t__, -1, 1> length_sigma =
         Eigen::Matrix<local_scalar_t__, -1, 1>::Constant(J, DUMMY_VAR__);
      for (int sym1__ = 1; sym1__ <= J; ++sym1__) {
        stan::model::assign(length_sigma, in__.read<local_scalar_t__>(),
          "assigning variable length_sigma", stan::model::index_uni(sym1__));
      }
      out__.write_free_lb(0, length_sigma);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
    } // transform_inits_impl() 
    
  inline void get_param_names(std::vector<std::string>& names__) const {
    
    names__ = std::vector<std::string>{"width_mu_0", "length_mu_0",
      "width_sigma_0", "length_sigma_0", "width_mu", "length_mu",
      "width_sigma", "length_sigma"};
    
    } // get_param_names() 
    
  inline void get_dims(std::vector<std::vector<size_t>>& dimss__) const {
    
    dimss__ = std::vector<std::vector<size_t>>{std::vector<size_t>{},
      std::vector<size_t>{}, std::vector<size_t>{}, std::vector<size_t>{
      }, std::vector<size_t>{static_cast<size_t>(J)},
      std::vector<size_t>{static_cast<size_t>(J)},
      std::vector<size_t>{static_cast<size_t>(J)},
      std::vector<size_t>{static_cast<size_t>(J)}};
    
    } // get_dims() 
    
  inline void constrained_param_names(
                                      std::vector<std::string>& param_names__,
                                      bool emit_transformed_parameters__ = true,
                                      bool emit_generated_quantities__ = true) const
    final {
    
    param_names__.emplace_back(std::string() + "width_mu_0");
    param_names__.emplace_back(std::string() + "length_mu_0");
    param_names__.emplace_back(std::string() + "width_sigma_0");
    param_names__.emplace_back(std::string() + "length_sigma_0");
    for (int sym1__ = 1; sym1__ <= J; ++sym1__) {
      {
        param_names__.emplace_back(std::string() + "width_mu" + '.' + std::to_string(sym1__));
      } 
    }
    for (int sym1__ = 1; sym1__ <= J; ++sym1__) {
      {
        param_names__.emplace_back(std::string() + "length_mu" + '.' + std::to_string(sym1__));
      } 
    }
    for (int sym1__ = 1; sym1__ <= J; ++sym1__) {
      {
        param_names__.emplace_back(std::string() + "width_sigma" + '.' + std::to_string(sym1__));
      } 
    }
    for (int sym1__ = 1; sym1__ <= J; ++sym1__) {
      {
        param_names__.emplace_back(std::string() + "length_sigma" + '.' + std::to_string(sym1__));
      } 
    }
    if (emit_transformed_parameters__) {
      
    }
    
    if (emit_generated_quantities__) {
      
    }
    
    } // constrained_param_names() 
    
  inline void unconstrained_param_names(
                                        std::vector<std::string>& param_names__,
                                        bool emit_transformed_parameters__ = true,
                                        bool emit_generated_quantities__ = true) const
    final {
    
    param_names__.emplace_back(std::string() + "width_mu_0");
    param_names__.emplace_back(std::string() + "length_mu_0");
    param_names__.emplace_back(std::string() + "width_sigma_0");
    param_names__.emplace_back(std::string() + "length_sigma_0");
    for (int sym1__ = 1; sym1__ <= J; ++sym1__) {
      {
        param_names__.emplace_back(std::string() + "width_mu" + '.' + std::to_string(sym1__));
      } 
    }
    for (int sym1__ = 1; sym1__ <= J; ++sym1__) {
      {
        param_names__.emplace_back(std::string() + "length_mu" + '.' + std::to_string(sym1__));
      } 
    }
    for (int sym1__ = 1; sym1__ <= J; ++sym1__) {
      {
        param_names__.emplace_back(std::string() + "width_sigma" + '.' + std::to_string(sym1__));
      } 
    }
    for (int sym1__ = 1; sym1__ <= J; ++sym1__) {
      {
        param_names__.emplace_back(std::string() + "length_sigma" + '.' + std::to_string(sym1__));
      } 
    }
    if (emit_transformed_parameters__) {
      
    }
    
    if (emit_generated_quantities__) {
      
    }
    
    } // unconstrained_param_names() 
    
  inline std::string get_constrained_sizedtypes() const {
    
    return std::string("[{\"name\":\"width_mu_0\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"length_mu_0\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"width_sigma_0\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"length_sigma_0\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"width_mu\",\"type\":{\"name\":\"vector\",\"length\":" + std::to_string(J) + "},\"block\":\"parameters\"},{\"name\":\"length_mu\",\"type\":{\"name\":\"vector\",\"length\":" + std::to_string(J) + "},\"block\":\"parameters\"},{\"name\":\"width_sigma\",\"type\":{\"name\":\"vector\",\"length\":" + std::to_string(J) + "},\"block\":\"parameters\"},{\"name\":\"length_sigma\",\"type\":{\"name\":\"vector\",\"length\":" + std::to_string(J) + "},\"block\":\"parameters\"}]");
    
    } // get_constrained_sizedtypes() 
    
  inline std::string get_unconstrained_sizedtypes() const {
    
    return std::string("[{\"name\":\"width_mu_0\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"length_mu_0\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"width_sigma_0\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"length_sigma_0\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"width_mu\",\"type\":{\"name\":\"vector\",\"length\":" + std::to_string(J) + "},\"block\":\"parameters\"},{\"name\":\"length_mu\",\"type\":{\"name\":\"vector\",\"length\":" + std::to_string(J) + "},\"block\":\"parameters\"},{\"name\":\"width_sigma\",\"type\":{\"name\":\"vector\",\"length\":" + std::to_string(J) + "},\"block\":\"parameters\"},{\"name\":\"length_sigma\",\"type\":{\"name\":\"vector\",\"length\":" + std::to_string(J) + "},\"block\":\"parameters\"}]");
    
    } // get_unconstrained_sizedtypes() 
    
  
    // Begin method overload boilerplate
    template <typename RNG>
    inline void write_array(RNG& base_rng,
                            Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                            Eigen::Matrix<double,Eigen::Dynamic,1>& vars,
                            const bool emit_transformed_parameters = true,
                            const bool emit_generated_quantities = true,
                            std::ostream* pstream = nullptr) const {
      const size_t num_params__ = 
  (((((((1 + 1) + 1) + 1) + J) + J) + J) + J);
      const size_t num_transformed = emit_transformed_parameters * 0;
      const size_t num_gen_quantities = emit_generated_quantities * 0;
      const size_t num_to_write = num_params__ + num_transformed +
        num_gen_quantities;
      std::vector<int> params_i;
      vars = Eigen::Matrix<double, Eigen::Dynamic, 1>::Constant(num_to_write,
        std::numeric_limits<double>::quiet_NaN());
      write_array_impl(base_rng, params_r, params_i, vars,
        emit_transformed_parameters, emit_generated_quantities, pstream);
    }

    template <typename RNG>
    inline void write_array(RNG& base_rng, std::vector<double>& params_r,
                            std::vector<int>& params_i,
                            std::vector<double>& vars,
                            bool emit_transformed_parameters = true,
                            bool emit_generated_quantities = true,
                            std::ostream* pstream = nullptr) const {
      const size_t num_params__ = 
  (((((((1 + 1) + 1) + 1) + J) + J) + J) + J);
      const size_t num_transformed = emit_transformed_parameters * 0;
      const size_t num_gen_quantities = emit_generated_quantities * 0;
      const size_t num_to_write = num_params__ + num_transformed +
        num_gen_quantities;
      vars = std::vector<double>(num_to_write,
        std::numeric_limits<double>::quiet_NaN());
      write_array_impl(base_rng, params_r, params_i, vars,
        emit_transformed_parameters, emit_generated_quantities, pstream);
    }

    template <bool propto__, bool jacobian__, typename T_>
    inline T_ log_prob(Eigen::Matrix<T_,Eigen::Dynamic,1>& params_r,
                       std::ostream* pstream = nullptr) const {
      Eigen::Matrix<int, -1, 1> params_i;
      return log_prob_impl<propto__, jacobian__>(params_r, params_i, pstream);
    }

    template <bool propto__, bool jacobian__, typename T__>
    inline T__ log_prob(std::vector<T__>& params_r,
                        std::vector<int>& params_i,
                        std::ostream* pstream = nullptr) const {
      return log_prob_impl<propto__, jacobian__>(params_r, params_i, pstream);
    }


    inline void transform_inits(const stan::io::var_context& context,
                         Eigen::Matrix<double, Eigen::Dynamic, 1>& params_r,
                         std::ostream* pstream = nullptr) const final {
      std::vector<double> params_r_vec(params_r.size());
      std::vector<int> params_i;
      transform_inits(context, params_i, params_r_vec, pstream);
      params_r = Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,1>>(
        params_r_vec.data(), params_r_vec.size());
    }

  inline void transform_inits(const stan::io::var_context& context,
                              std::vector<int>& params_i,
                              std::vector<double>& vars,
                              std::ostream* pstream__ = nullptr) const {
     constexpr std::array<const char*, 8> names__{"width_mu_0",
      "length_mu_0", "width_sigma_0", "length_sigma_0", "width_mu",
      "length_mu", "width_sigma", "length_sigma"};
      const std::array<Eigen::Index, 8> constrain_param_sizes__{1, 1, 
       1, 1, J, J, J, J};
      const auto num_constrained_params__ = std::accumulate(
        constrain_param_sizes__.begin(), constrain_param_sizes__.end(), 0);
    
     std::vector<double> params_r_flat__(num_constrained_params__);
     Eigen::Index size_iter__ = 0;
     Eigen::Index flat_iter__ = 0;
     for (auto&& param_name__ : names__) {
       const auto param_vec__ = context.vals_r(param_name__);
       for (Eigen::Index i = 0; i < constrain_param_sizes__[size_iter__]; ++i) {
         params_r_flat__[flat_iter__] = param_vec__[i];
         ++flat_iter__;
       }
       ++size_iter__;
     }
     vars.resize(num_params_r__);
     transform_inits_impl(params_r_flat__, params_i, vars, pstream__);
    } // transform_inits() 
    
};
}

using stan_model = hierarchical_size_model_namespace::hierarchical_size_model;

#ifndef USING_R

// Boilerplate
stan::model::model_base& new_model(
        stan::io::var_context& data_context,
        unsigned int seed,
        std::ostream* msg_stream) {
  stan_model* m = new stan_model(data_context, seed, msg_stream);
  return *m;
}

stan::math::profile_map& get_stan_profile_data() {
  return hierarchical_size_model_namespace::profiles__;
}

#endif

