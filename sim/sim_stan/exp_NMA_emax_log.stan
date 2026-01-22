// AQSS 

functions {
 /* compute correlated group-level effects
  * in the presence of a within-group covariance matrix
  * Args:
  *   z: matrix of unscaled group-level effects
  *   SD: vector of standard deviation parameters
  *   L: cholesky factor correlation matrix
  *   Lcov: cholesky factor of within-group correlation matrix
  * Returns:
  *   matrix of scaled group-level effects
  */
  matrix scale_r_cor_cov(matrix z, vector SD, matrix L, matrix Lcov) {
    vector[num_elements(z)] z_flat = to_vector(z);
    vector[num_elements(z)] r = rep_vector(0, num_elements(z));
    matrix[rows(L), cols(L)] LC = diag_pre_multiply(SD, L);
    int rows_z = rows(z);
    int rows_L = rows(L);
    // kronecker product of cholesky factors times a vector
    for (icov in 1:rows(Lcov)) {
      for (jcov in 1:icov) {
        if (Lcov[icov, jcov] > 1e-10) {
          // avoid calculating products between unrelated individuals
          for (i in 1:rows_L) {
            for (j in 1:i) {
              // incremented element of the output vector
              int k = (rows_L * (icov - 1)) + i;
              // applied element of the input vector
              int l = (rows_L * (jcov - 1)) + j;
              r[k] = r[k] + Lcov[icov, jcov] * LC[i, j] * z_flat[l];
            }
          }
        }
      }
    }
    // r is returned in another dimension order than z
    return to_matrix(r, cols(z), rows(z), 0);
  }
 // ------------------- Additional stan functions    
  /* multi-normal log-PDF for fixed correlation matrices
   * assuming homogoneous variances
   * Args:
   *   y: response vector
   *   mu: mean parameter vector
   *   sigma: residual standard deviation
   *   chol_cor: cholesky factor of the correlation matrix
   * Returns:
   *   sum of the log-PDF values of all observations
   */
  real normal_fcor_hom_lpdf(vector y, vector mu, real sigma, data matrix chol_cor) {
    return multi_normal_cholesky_lpdf(y | mu, sigma * chol_cor);
  }
  /* multi-normal log-PDF for fixed correlation matrices
   * assuming heterogenous variances
   * Args:
   *   y: response vector
   *   mu: mean parameter vector
   *   sigma: residual standard deviation vector
   *   chol_cor: cholesky factor of the correlation matrix
   * Returns:
   *   sum of the log-PDF values of all observations
   */
  real normal_fcor_het_lpdf(vector y, vector mu, vector sigma, data matrix chol_cor) {
    return multi_normal_cholesky_lpdf(y | mu, diag_pre_multiply(sigma, chol_cor));
  }
  /* multi-student-t log-PDF for fixed correlation matrices
   * assuming homogoneous variances
   * Args:
   *   y: response vector
   *   nu: degrees of freedom parameter
   *   mu: mean parameter vector
   *   sigma: scale parameter
   *   chol_cor: cholesky factor of the correlation matrix
   * Returns:
   *   sum of the log-PDF values of all observations
   */
  real student_t_fcor_hom_lpdf(vector y, real nu, vector mu, real sigma,
                               data matrix chol_cor) {
    int N = rows(chol_cor);
    matrix[N, N] Cov = multiply_lower_tri_self_transpose(sigma * chol_cor);
    return multi_student_t_lpdf(y | nu, mu, Cov);
  }
  /* multi-student-t log-PDF for fixed correlation matrices
   * assuming heterogenous variances
   * Args:
   *   y: response vector
   *   nu: degrees of freedom parameter
   *   mu: mean parameter vector
   *   sigma: scale parameter vector
   *   chol_cor: cholesky factor of the correlation matrix
   * Returns:
   *   sum of the log-PDF values of all observations
   */
  real student_t_fcor_het_lpdf(vector y, real nu, vector mu, vector sigma,
                               data matrix chol_cor) {
    int N = rows(chol_cor);
    matrix[N, N] Cov = diag_pre_multiply(sigma, chol_cor);
    Cov = multiply_lower_tri_self_transpose(Cov);
    return multi_student_t_lpdf(y | nu, mu, Cov);
  }
 // ------------------- end additional functions  
} // end function block

data {
  int<lower=1> N;  // total number of observations
  vector[N] Y;  // response variable
  vector<lower=0>[N] se;  // known sampling error
  matrix[N, N] Mfcor;  // known residual covariance matrix
  int<lower=1> K_MUemaxT;  // number of population-level effects
  matrix[N, K_MUemaxT] X_MUemaxT;  // population-level design matrix
  int<lower=1> K_emaxA;  // number of population-level effects
  matrix[N, K_emaxA] X_emaxA;  // population-level design matrix
  int<lower=1> K_emaxB;  // number of population-level effects
  matrix[N, K_emaxB] X_emaxB;  // population-level design matrix
  int<lower=1> K_emaxC;  // number of population-level effects
  matrix[N, K_emaxC] X_emaxC;  // population-level design matrix
  int<lower=1> K_ed50A;  // number of population-level effects
  matrix[N, K_ed50A] X_ed50A;  // population-level design matrix
  int<lower=1> K_ed50B;  // number of population-level effects
  matrix[N, K_ed50B] X_ed50B;  // population-level design matrix
  int<lower=1> K_ed50C;  // number of population-level effects
  matrix[N, K_ed50C] X_ed50C;  // population-level design matrix
  int<lower=1> K_MUkT;  // number of population-level effects
  matrix[N, K_MUkT] X_MUkT;  // population-level design matrix
  int<lower=1> K_sA;  // number of population-level effects
  matrix[N, K_sA] X_sA;  // population-level design matrix
  int<lower=1> K_sB;  // number of population-level effects
  matrix[N, K_sB] X_sB;  // population-level design matrix
  int<lower=1> K_sC;  // number of population-level effects
  matrix[N, K_sC] X_sC;  // population-level design matrix
  // covariates for non-linear functions
  vector[N] C_1;
  vector[N] C_2;
  vector[N] C_3;
  vector[N] C_4;
  vector[N] C_5;
  // data for group-level effects of ID 1
  int<lower=1> N_1;  // number of grouping levels
  int<lower=1> M_1;  // number of coefficients per level
  array[N] int<lower=1> J_1;  // grouping indicator per observation
  matrix[N_1, N_1] Lcov_1;  // cholesky factor of known covariance matrix
  // group-level predictor values
  vector[N] Z_1_MUemaxT_1;
  vector[N] Z_1_MUkT_2;
  int<lower=1> NC_1;  // number of group-level correlations
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
  matrix[N, N] Lfcor = cholesky_decompose(Mfcor); // Add cholesky decomposition
}
parameters {
  vector[K_MUemaxT] b_MUemaxT;  // regression coefficients
  vector[K_emaxA] b_emaxA;  // regression coefficients
  vector[K_emaxB] b_emaxB;  // regression coefficients
  vector[K_emaxC] b_emaxC;  // regression coefficients
  vector[K_ed50A] b_ed50A;  // regression coefficients
  vector[K_ed50B] b_ed50B;  // regression coefficients
  vector[K_ed50C] b_ed50C;  // regression coefficients
  vector[K_MUkT] b_MUkT;  // regression coefficients
  vector[K_sA] b_sA;  // regression coefficients
  vector[K_sB] b_sB;  // regression coefficients
  vector[K_sC] b_sC;  // regression coefficients
  vector<lower=0>[M_1] sd_1;  // group-level standard deviations
  matrix[M_1, N_1] z_1;  // standardized group-level effects
  cholesky_factor_corr[M_1] L_1;  // cholesky factor of correlation matrix
}
transformed parameters {
  real sigma = 1;  // change setting of dispersion parameter to 1
  matrix[N_1, M_1] r_1;  // actual group-level effects
  // using vectors speeds up indexing in loops
  vector[N_1] r_1_MUemaxT_1;
  vector[N_1] r_1_MUkT_2;
  real lprior = 0;  // prior contributions to the log posterior
  // compute actual group-level effects
  r_1 = scale_r_cor_cov(z_1, sd_1, L_1, Lcov_1);
  r_1_MUemaxT_1 = r_1[, 1];
  r_1_MUkT_2 = r_1[, 2];
  lprior += normal_lpdf(b_MUemaxT[1] | -2.5, 10);
  lprior += normal_lpdf(b_MUemaxT[2] | -2.5, 10);
  lprior += normal_lpdf(b_MUemaxT[3] | -2.5, 10);
  lprior += normal_lpdf(b_MUemaxT[4] | -2.5, 10);
  lprior += normal_lpdf(b_MUemaxT[5] | -2.5, 10);
  lprior += normal_lpdf(b_MUemaxT[6] | -2.5, 10);
  lprior += normal_lpdf(b_MUemaxT[7] | -2.5, 10);
  lprior += normal_lpdf(b_MUemaxT[8] | -2.5, 10);
  lprior += normal_lpdf(b_MUemaxT[9] | -2.5, 10);
  lprior += normal_lpdf(b_MUemaxT[10] | -2.5, 10);
  lprior += normal_lpdf(b_emaxA | -10, 10);
  lprior += normal_lpdf(b_emaxB | -15, 10);
  lprior += normal_lpdf(b_emaxC | -20, 10);
  lprior += normal_lpdf(b_ed50A | -2, 10);
  lprior += normal_lpdf(b_ed50B | -1.5, 10);
  lprior += normal_lpdf(b_ed50C | -0.75, 10);
  lprior += normal_lpdf(b_MUkT[1] | -2.5, 10);
  lprior += normal_lpdf(b_MUkT[2] | -2.5, 10);
  lprior += normal_lpdf(b_MUkT[3] | -2.5, 10);
  lprior += normal_lpdf(b_MUkT[4] | -2.5, 10);
  lprior += normal_lpdf(b_MUkT[5] | -2.5, 10);
  lprior += normal_lpdf(b_MUkT[6] | -2.5, 10);
  lprior += normal_lpdf(b_MUkT[7] | -2.5, 10);
  lprior += normal_lpdf(b_MUkT[8] | -2.5, 10);
  lprior += normal_lpdf(b_MUkT[9] | -2.5, 10);
  lprior += normal_lpdf(b_MUkT[10] | -2.5, 10);
  lprior += normal_lpdf(b_sA | 0, 10);
  lprior += normal_lpdf(b_sB | 0, 10);
  lprior += normal_lpdf(b_sC | 0, 10);
  lprior += normal_lpdf(sd_1 | 0, 5)
    - 2 * normal_lccdf(0 | 0, 5);
  lprior += lkj_corr_cholesky_lpdf(L_1 | 1);
}
model {
  // likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N] nlp_MUemaxT = rep_vector(0.0, N);
    // initialize linear predictor term
    vector[N] nlp_emaxA = rep_vector(0.0, N);
    // initialize linear predictor term
    vector[N] nlp_emaxB = rep_vector(0.0, N);
    // initialize linear predictor term
    vector[N] nlp_emaxC = rep_vector(0.0, N);
    // initialize linear predictor term
    vector[N] nlp_ed50A = rep_vector(0.0, N);
    // initialize linear predictor term
    vector[N] nlp_ed50B = rep_vector(0.0, N);
    // initialize linear predictor term
    vector[N] nlp_ed50C = rep_vector(0.0, N);
    // initialize linear predictor term
    vector[N] nlp_MUkT = rep_vector(0.0, N);
    // initialize linear predictor term
    vector[N] nlp_sA = rep_vector(0.0, N);
    // initialize linear predictor term
    vector[N] nlp_sB = rep_vector(0.0, N);
    // initialize linear predictor term
    vector[N] nlp_sC = rep_vector(0.0, N);
    // initialize non-linear predictor term
    vector[N] mu;
    nlp_MUemaxT += X_MUemaxT * b_MUemaxT;
    nlp_emaxA += X_emaxA * b_emaxA;
    nlp_emaxB += X_emaxB * b_emaxB;
    nlp_emaxC += X_emaxC * b_emaxC;
    nlp_ed50A += X_ed50A * b_ed50A;
    nlp_ed50B += X_ed50B * b_ed50B;
    nlp_ed50C += X_ed50C * b_ed50C;
    nlp_MUkT += X_MUkT * b_MUkT;
    nlp_sA += X_sA * b_sA;
    nlp_sB += X_sB * b_sB;
    nlp_sC += X_sC * b_sC;
    for (n in 1:N) {
      // add more terms to the linear predictor
      nlp_MUemaxT[n] += r_1_MUemaxT_1[J_1[n]] * Z_1_MUemaxT_1[n];
    }
    for (n in 1:N) {
      // add more terms to the linear predictor
      nlp_MUkT[n] += r_1_MUkT_2[J_1[n]] * Z_1_MUkT_2[n];
    }
    for (n in 1:N) {
      // compute non-linear predictor values
      mu[n] = ((nlp_MUemaxT[n] + (nlp_emaxA[n] * C_1[n] + nlp_emaxB[n] * C_2[n] + nlp_emaxC[n] * C_3[n]) * C_4[n] / (exp(nlp_ed50A[n] * C_1[n] + nlp_ed50B[n] * C_2[n] + nlp_ed50C[n] * C_3[n]) + C_4[n])) * (1 - exp( - exp(nlp_MUkT[n] + nlp_sA[n] * C_1[n] * log(C_4[n] + 1) + nlp_sB[n] * C_2[n] * log(C_4[n] + 1) + nlp_sC[n] * C_3[n] * log(C_4[n] + 1)) * C_5[n])));
    }
    target += normal_fcor_hom_lpdf(Y | mu, sigma, Lfcor); // Change log-PDF
  }
  // priors including constants
  target += lprior;
  target += std_normal_lpdf(to_vector(z_1));
}
generated quantities {
  // compute group-level correlations
  corr_matrix[M_1] Cor_1 = multiply_lower_tri_self_transpose(L_1);
  vector<lower=-1,upper=1>[NC_1] cor_1;
  // extract upper diagonal of correlation matrix
  for (k in 1:M_1) {
    for (j in 1:(k - 1)) {
      cor_1[choose(k - 1, 2) + j] = Cor_1[j, k];
    }
  }
}


