#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/BLAS.h> 
#include <R_ext/Lapack.h>
#include <R_ext/Rdynload.h>

// Auxiliar functions for fllp_qmc

static inline double norm_inv_cdf(double p) {
  return qnorm(p, 0.0, 1.0, 1, 0);
}

static inline double norm_cdf(double x) {
  return pnorm(x, 0.0, 1.0, 1, 0);
}

// Auxiliar functions for rtmvnorm_gibbs

// truncated normal univariate
static double rtnorm_scalar(double mean, double sd, double lower, double upper) {
  double Z1 = (lower - mean) / sd;
  double Z2 = (upper - mean) / sd;

  double pl = pnorm(Z1, 0.0, 1.0, 1, 0);
  double pu = pnorm(Z2, 0.0, 1.0, 1, 0);

  double u = unif_rand();
  double p = (pu - pl) * u + pl;

  if (p < 1e-15) p = 1e-15;
  if (p > 1 - 1e-15) p = 1 - 1e-15;

  return mean + sd * qnorm(p, 0.0, 1.0, 1, 0);
}

// vector substraction j (dim n-1)
static void vector_without_j(const double *vec, double *out, int n, int j) {
  int pos = 0;
  for (int i = 0; i < n; i++) {
    if (i != j) out[pos++] = vec[i];
  }
}

// sub-matriz n x n without row and column j (dim n-1 x n-1)
static void submatrix_without_j(const double *mat, double *out, int n, int j) {
  int col_out = 0;
  for (int col = 0; col < n; col++) {
    if (col == j) continue;
    int row_out = 0;
    for (int row = 0; row < n; row++) {
      if (row == j) continue;
      out[col_out * (n - 1) + row_out] = mat[col * n + row];
      row_out++;
    }
    col_out++;
  }
}

// truncated normal multivariate
SEXP rtmvnorm_gibbs(SEXP mu_, SEXP lb_, SEXP ub_, SEXP Sigma_, SEXP V_,
                    SEXP n_, SEXP ntraits_,
                    SEXP n_iter_, SEXP burn_in_, SEXP thin_) {

  PROTECT(mu_ = coerceVector(mu_, REALSXP));
  PROTECT(lb_ = coerceVector(lb_, REALSXP));
  PROTECT(ub_ = coerceVector(ub_, REALSXP));
  PROTECT(Sigma_ = coerceVector(Sigma_, REALSXP));
  PROTECT(V_ = coerceVector(V_, REALSXP));
  PROTECT(n_ = coerceVector(n_, INTSXP));
  PROTECT(ntraits_ = coerceVector(ntraits_, INTSXP));
  PROTECT(n_iter_ = coerceVector(n_iter_, INTSXP));
  PROTECT(burn_in_ = coerceVector(burn_in_, INTSXP));
  PROTECT(thin_ = coerceVector(thin_, INTSXP));

  int n = INTEGER(n_)[0];
  int ntraits = INTEGER(ntraits_)[0];
  int n_iter = INTEGER(n_iter_)[0];
  int burn_in = INTEGER(burn_in_)[0];
  int thin = INTEGER(thin_)[0];

  double *mu = REAL(mu_);
  double *lb = REAL(lb_);
  double *ub = REAL(ub_);
  double *Sigma = REAL(Sigma_);
  double *V = REAL(V_);

  int inc = 1;
  char trans = 'T';
  double alpha = 1.0, beta = 0.0;

  SEXP out_;
  PROTECT(out_ = allocMatrix(REALSXP, n, ntraits));
  double *out = REAL(out_);

  GetRNGstate();

  int dim_sub = ntraits - 1;
  double *V_j_minus_j = (double *) R_alloc(dim_sub, sizeof(double));
  double *Sigma_j_minus_j = (double *) R_alloc(dim_sub, sizeof(double));
  double *V_minus_j_minus_j = (double *) R_alloc(dim_sub * dim_sub, sizeof(double));
  double *SigmaInv_jj = (double *) R_alloc(dim_sub * dim_sub, sizeof(double));
  double *temp_vec = (double *) R_alloc(dim_sub, sizeof(double));
  double *Aj = (double *) R_alloc(dim_sub, sizeof(double));

  double *A_list = (double *) R_alloc(ntraits * dim_sub, sizeof(double));
  double *cond_var_list = (double *) R_alloc(ntraits, sizeof(double));

  // Pre-compute A_list and cond_var_list
  for (int j = 0; j < ntraits; j++) {
    double *V_j = V + (long long)j * ntraits;
    double V_jj = V_j[j];

    vector_without_j(V_j, V_j_minus_j, ntraits, j);
    submatrix_without_j(V, V_minus_j_minus_j, ntraits, j);

    // SigmaInv_jj = V_minus_j_minus_j - (V_j_minus_j %*% t(V_j_minus_j)) / V_jj
    memcpy(SigmaInv_jj, V_minus_j_minus_j, dim_sub * dim_sub * sizeof(double));
    double minus_inv_alpha = -1.0 / V_jj;
    F77_NAME(dger)(&dim_sub, &dim_sub, &minus_inv_alpha,
                 V_j_minus_j, &inc,
                 V_j_minus_j, &inc,
                 SigmaInv_jj, &dim_sub);

    double *Sigma_j = Sigma + (long long)j * ntraits;
    vector_without_j(Sigma_j, Sigma_j_minus_j, ntraits, j);

    F77_NAME(dgemv)(&trans, &dim_sub, &dim_sub, &alpha, SigmaInv_jj, &dim_sub, Sigma_j_minus_j, &inc, &beta, Aj, &inc FCONE);
    for (int k = 0; k < dim_sub; k++) {
      A_list[j * dim_sub + k] = Aj[k];
    }

    cond_var_list[j] = Sigma[j * ntraits + j] - F77_NAME(ddot)(&dim_sub, Aj, &inc, Sigma_j_minus_j, &inc);
  }

  double *mu_k = (double *) R_alloc(ntraits, sizeof(double));
  double *lb_k = (double *) R_alloc(ntraits, sizeof(double));
  double *ub_k = (double *) R_alloc(ntraits, sizeof(double));
  double *x = (double *) R_alloc(ntraits, sizeof(double));
  double *mean_samples = (double *) R_alloc(ntraits, sizeof(double));
  double *x_minus_j = (double *) R_alloc(dim_sub, sizeof(double));
  double *mu_minus_j = (double *) R_alloc(dim_sub, sizeof(double));

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < ntraits; j++) {
      mu_k[j] = mu[j * n + i];
      lb_k[j] = lb[j * n + i];
      ub_k[j] = ub[j * n + i];
      x[j] = mu_k[j];  // inicialization
      mean_samples[j] = 0.0;
    }
    int sample_count = 1;

    for (int iter = 0; iter < n_iter; iter++) {
      for (int j = 0; j < ntraits; j++) {
        int pos = 0;
        for (int idx = 0; idx < ntraits; idx++) {
          if (idx != j) {
            x_minus_j[pos] = x[idx];
            mu_minus_j[pos] = mu_k[idx];
            pos++;
          }
        }
        for (int k = 0; k < dim_sub; k++) temp_vec[k] = x_minus_j[k] - mu_minus_j[k];
        double cond_mean = mu_k[j] + F77_NAME(ddot)(&dim_sub, A_list + j * dim_sub, &inc, temp_vec, &inc);
        double cond_sd = sqrt(cond_var_list[j]);
        x[j] = rtnorm_scalar(cond_mean, cond_sd, lb_k[j], ub_k[j]);
      }
      if (iter >= burn_in && (iter % thin == 0)) {
        for (int j = 0; j < ntraits; j++) {
          mean_samples[j] = (x[j] + mean_samples[j] * (sample_count - 1)) / sample_count;
        }
        sample_count++;
      }
    }
    for (int j = 0; j < ntraits; j++) {
      out[j * n + i] = mean_samples[j];
    }
  }

  PutRNGstate();

  UNPROTECT(11);
  return out_;
}

SEXP rmvn(SEXP n_, SEXP ntraits_, SEXP mean_, SEXP Sigma_) {
  PROTECT(n_ = coerceVector(n_, INTSXP));
  PROTECT(ntraits_ = coerceVector(ntraits_, INTSXP));
  PROTECT(mean_ = coerceVector(mean_, REALSXP));
  PROTECT(Sigma_ = coerceVector(Sigma_, REALSXP));

  int n = INTEGER(n_)[0];
  int ntraits = INTEGER(ntraits_)[0];
  double *mean = REAL(mean_);
  double *Sigma = REAL(Sigma_);

  SEXP out_;
  PROTECT(out_ = allocMatrix(REALSXP, n, ntraits));
  double *out = REAL(out_);

  double *L = (double *) R_alloc(ntraits * ntraits, sizeof(double));
  memcpy(L, Sigma, ntraits * ntraits * sizeof(double));

  int info;
  F77_NAME(dpotrf)("L", &ntraits, L, &ntraits, &info FCONE);  // Cholesky
  if (info != 0) error("Cholesky factorization failed: dpotrf returned %d", info);

  double *z = (double *) R_alloc(ntraits, sizeof(double));
  double *x = (double *) R_alloc(ntraits, sizeof(double));
  char trans = 'N';
  double alpha = 1.0, beta = 0.0;
  int inc = 1;

  GetRNGstate();

  for (int i = 0; i < n; i++) {

    for (int j = 0; j < ntraits; j++) {
      z[j] = norm_rand();  // z ~ N(0, 1)
    }

    // x = L * z

    F77_NAME(dgemv)(&trans, &ntraits, &ntraits, &alpha,
                    L, &ntraits, z, &inc, &beta, x, &inc FCONE);

    // out[i, ] = mean[i, ] + x
    for (int j = 0; j < ntraits; j++) {
      out[i + j * n] = mean[i + j * n] + x[j]; 
    }
  }

  PutRNGstate();
  UNPROTECT(5);
  return out_;
}

SEXP fBj_DLN_G_mtme(SEXP Bv, SEXP SigmaBInv, SEXP SigmaEInv,
           SEXP rv, SEXP X, SEXP x2, SEXP n, SEXP p, SEXP ntraits) {

  PROTECT(Bv = coerceVector(Bv, REALSXP));                double* pBv = REAL(Bv);
  PROTECT(SigmaBInv = coerceVector(SigmaBInv, REALSXP));  double* pSigmaBInv = REAL(SigmaBInv);
  PROTECT(SigmaEInv = coerceVector(SigmaEInv, REALSXP));  double* pSigmaEInv = REAL(SigmaEInv);
  PROTECT(rv = coerceVector(rv, REALSXP));                double* prv = REAL(rv);
  PROTECT(X = coerceVector(X, REALSXP));                  double* pX = REAL(X);
  PROTECT(x2 = coerceVector(x2, REALSXP));                double* px2 = REAL(x2);
  PROTECT(n = coerceVector(n, INTSXP));                   int N = INTEGER(n)[0];
  PROTECT(p = coerceVector(p, INTSXP));                   int P = INTEGER(p)[0];
  PROTECT(ntraits = coerceVector(ntraits, INTSXP));       int T = INTEGER(ntraits)[0];

  SEXP result;
  PROTECT(result = allocVector(VECSXP, 2));

  int inc = 1;
  int info;

  char transN = 'N', transT = 'T', uploL = 'L';
  double one = 1.0, zero = 0.0;

  //------------------------------------------------------------------
  // Detectar si SigmaBInv y SigmaEInv son diagonales
  //------------------------------------------------------------------

  bool diagonalCase = (T == 1);

  if (T > 1) {
    diagonalCase = true;

    for (int j = 0; j < T && diagonalCase; ++j) {
      for (int i = 0; i < T; ++i) {

        if (i == j)
          continue;

        if (pSigmaBInv[i + j * T] != 0.0 ||
            pSigmaEInv[i + j * T] != 0.0) {
          diagonalCase = false;
          break;
        }
      }
    }
  }

  //------------------------------------------------------------------
  // Memoria auxiliar
  //------------------------------------------------------------------

  double *Xk;
  double *Xt_rv = (double*) R_alloc(T, sizeof(double));
  double *Bkt   = (double*) R_alloc(T, sizeof(double));

  double *SigmaBk = nullptr;
  double *tmp     = nullptr;
  double *z       = nullptr;

  if (!diagonalCase) {
    SigmaBk = (double*) R_alloc(T * T, sizeof(double));
    tmp     = (double*) R_alloc(T, sizeof(double));
    z       = (double*) R_alloc(T, sizeof(double));
  }

  GetRNGstate();

  for (int k = 0; k < P; ++k) {

    Xk = pX + (long long)k * N;

    //----------------------------------------------------------------
    // CASO ESCALAR
    //----------------------------------------------------------------

    if (T == 1) {

      double Bk = pBv[k];

      double rhs =
        F77_CALL(ddot)(
          &N,
          Xk, &inc,
          prv, &inc);

      rhs += px2[k] * Bk;

      double c =
        pSigmaBInv[0] +
        px2[k] * pSigmaEInv[0];

      double mu =
        pSigmaEInv[0] * rhs / c;

      pBv[k] =
        mu +
        norm_rand() / sqrt(c);

      double diff = Bk - pBv[k];

      F77_CALL(daxpy)(
        &N,
        &diff,
        Xk,
        &inc,
        prv,
        &inc);

      continue;
    }

    //----------------------------------------------------------------
    // CASO DIAGONAL
    //----------------------------------------------------------------

    if (diagonalCase) {

      for (int t = 0; t < T; ++t) {

        Bkt[t] = pBv[k + t * P];

        double rhs =
          F77_CALL(ddot)(
            &N,
            Xk, &inc,
            prv + t * N, &inc);

        rhs += px2[k] * Bkt[t];

        double c =
          pSigmaBInv[t + t * T] +
          px2[k] * pSigmaEInv[t + t * T];

        double mu =
          pSigmaEInv[t + t * T] *
          rhs / c;

        pBv[k + t * P] =
          mu +
          norm_rand() / sqrt(c);
      }

      for (int t = 0; t < T; ++t) {

        double diff =
          Bkt[t] -
          pBv[k + t * P];

        F77_CALL(daxpy)(
          &N,
          &diff,
          Xk,
          &inc,
          prv + t * N,
          &inc);
      }

      continue;
    }

    //----------------------------------------------------------------
    // CASO GENERAL
    //----------------------------------------------------------------

    // Xt_rv = t(Xk) %*% rv + x2[k] * Bk
    for (int t = 0; t < T; ++t) {

      Bkt[t] = pBv[k + t * P];

      Xt_rv[t] =
        F77_CALL(ddot)(
          &N,
          Xk,
          &inc,
          prv + t * N,
          &inc)
        + px2[k] * Bkt[t];
    }

    // SigmaBk = SigmaBInv + x2[k] * SigmaEInv
    for (int i = 0; i < T * T; ++i)
      SigmaBk[i] =
        pSigmaBInv[i] +
        px2[k] * pSigmaEInv[i];

    // Cholesky: SigmaBk = L L'
    F77_CALL(dpotrf)(
      &uploL,
      &T,
      SigmaBk,
      &T,
      &info FCONE);

    if (info != 0)
      error("Cholesky factorization failed (k=%d)", k);

    // tmp = SigmaEInv %*% Xt_rv
    F77_CALL(dgemv)(
      &transN,
      &T,
      &T,
      &one,
      pSigmaEInv,
      &T,
      Xt_rv,
      &inc,
      &zero,
      tmp,
      &inc FCONE);

    // tmp = solve(SigmaBk, tmp)
    F77_CALL(dpotrs)(
      &uploL,
      &T,
      &inc,
      SigmaBk,
      &T,
      tmp,
      &T,
      &info FCONE);

    if (info != 0)
      error("Solve linear system failed (k=%d)", k);

    // z ~ N(0, I)
    for (int i = 0; i < T; ++i)
      z[i] = norm_rand();

    // z <- L^{-T} z
    F77_CALL(dtrsv)(
      &uploL,
      &transT,
      "N",
      &T,
      SigmaBk,
      &T,
      z,
      &inc
      FCONE FCONE FCONE);

    for (int t = 0; t < T; ++t) {

      pBv[k + t * P] =
        tmp[t] + z[t];

      double diff =
        Bkt[t] -
        pBv[k + t * P];

      F77_CALL(daxpy)(
        &N,
        &diff,
        Xk,
        &inc,
        prv + t * N,
        &inc);
    }
  }

  PutRNGstate();

  SET_VECTOR_ELT(result, 0, Bv);
  SET_VECTOR_ELT(result, 1, rv);

  UNPROTECT(10);
  return result;
}

SEXP fBj_P_mtme(SEXP w, SEXP yr, SEXP Bv, SEXP SigmaBInv, SEXP rv, 
           SEXP X, SEXP X2, SEXP p, SEXP ntraits, SEXP n) {

  PROTECT(w           = coerceVector(w, REALSXP));           double *pw = REAL(w);
  PROTECT(yr          = coerceVector(yr, REALSXP));          double *pyr = REAL(yr);
  PROTECT(Bv          = coerceVector(Bv, REALSXP));          double *pBv = REAL(Bv);
  PROTECT(SigmaBInv   = coerceVector(SigmaBInv, REALSXP));   double *pSigmaBInv = REAL(SigmaBInv);
  PROTECT(rv          = coerceVector(rv, REALSXP));          double *prv = REAL(rv);
  PROTECT(X           = coerceVector(X, REALSXP));           double *pX = REAL(X);
  PROTECT(X2          = coerceVector(X2, REALSXP));          double *pX2 = REAL(X2);
  PROTECT(p           = coerceVector(p, INTSXP));            int P = INTEGER(p)[0];
  PROTECT(ntraits     = coerceVector(ntraits, INTSXP));      int T = INTEGER(ntraits)[0];
  PROTECT(n           = coerceVector(n, INTSXP));            int N = INTEGER(n)[0];

  int inc = 1;
  char transT = 'T', uploL = 'L';
  int info;
  double *Xk;
  double *X2k;

  double *omega_k = (double *) R_alloc(T, sizeof(double));
  double *SigmaBk = (double *) R_alloc(T * T, sizeof(double));
  double *score = (double *) R_alloc(T, sizeof(double));
  double *z = (double *) R_alloc(T, sizeof(double));
  double *res_t = (double *) R_alloc(N, sizeof(double));
  double* Bkt = (double*) R_alloc(T, sizeof(double));

  GetRNGstate();

  for (int k = 0; k < P; ++k) {
    Xk = pX + (long long)k * N;
    X2k = pX2 + (long long)k * N;

    // omega_k[t] = sum_i w[i + t*N] * X2[i + k*N] and SigmaBk = SigmaBInv + diag(omega_k)
    memcpy(SigmaBk, pSigmaBInv, T * T * sizeof(double));
    for (int t = 0; t < T; ++t) {
      Bkt[t] = -pBv[k + t * P];
      omega_k[t] = F77_CALL(ddot)(&N, pw + t * N, &inc, X2k, &inc); 
      SigmaBk[t + t * T] += omega_k[t];
    }

    // Cholesky factorization: SigmaBk = L L^T
    F77_CALL(dpotrf)(&uploL, &T, SigmaBk, &T, &info FCONE);
    if (info != 0) error("Cholesky failed in SigmaBk (k=%d)", k);

    // score[t] = sum_i Xk[i] * (yr[i + t*N] - w[i + t*N] * rv[i + t*N])
    for (int t = 0; t < T; ++t) {
      for (int i = 0; i < N; ++i) {
        res_t[i] = pyr[i + t * N] - pw[i + t * N] * (prv[i + t * N] + Xk[i]*Bkt[t]);
      }
      score[t] = F77_CALL(ddot)(&N, Xk, &inc, res_t, &inc);
    }

    // Solve SigmaBk * muBk = score
    F77_CALL(dpotrs)(&uploL, &T, &inc, SigmaBk, &T, score, &T, &info FCONE);
    if (info != 0) error("Solve linear system failed (k=%d)", k);
   
    // z ~ N(0, I)
    for (int t = 0; t < T; ++t)
      z[t] = norm_rand();

    // Solve: L^T y = z (i.e., z := L^{-T} z)
    F77_CALL(dtrsv)(&uploL, &transT, "N", &T, SigmaBk, &T, z, &inc FCONE FCONE FCONE);

    // rv += Xk * Bv.row(k)
    for (int t = 0; t < T; ++t) {
      // Bv.row(k) = muBk + z  
      pBv[k + t * P] = score[t] + z[t];
      double diff = pBv[k + t * P] + Bkt[t];
      F77_CALL(daxpy)(&N, &diff, Xk, &inc, prv + t * N, &inc);
    }
  }

  PutRNGstate();

  SEXP result = PROTECT(allocVector(VECSXP, 2));
  SET_VECTOR_ELT(result, 0, Bv);
  SET_VECTOR_ELT(result, 1, rv);
  UNPROTECT(11);
  return result;
}

SEXP fUi_mtme(SEXP w, SEXP yr, SEXP U, SEXP SigmaUInv, SEXP rv, 
         SEXP ntraits, SEXP n) {

  PROTECT(w           = coerceVector(w, REALSXP));           double *pw = REAL(w);
  PROTECT(yr          = coerceVector(yr, REALSXP));          double *pyr = REAL(yr);
  PROTECT(U           = coerceVector(U, REALSXP));           double *pU = REAL(U);
  PROTECT(SigmaUInv   = coerceVector(SigmaUInv, REALSXP));   double *pSigmaUInv = REAL(SigmaUInv);
  PROTECT(rv          = coerceVector(rv, REALSXP));          double *prv = REAL(rv);
  PROTECT(ntraits     = coerceVector(ntraits, INTSXP));      int T = INTEGER(ntraits)[0];
  PROTECT(n           = coerceVector(n, INTSXP));            int N = INTEGER(n)[0];

  int inc = 1;
  char transT = 'T', uploL = 'L';
  int info;

  double *omega_k = (double *) R_alloc(T, sizeof(double));
  double *SigmaUk = (double *) R_alloc(T * T, sizeof(double));
  double *score = (double *) R_alloc(T, sizeof(double));
  double *z = (double *) R_alloc(T, sizeof(double));

  GetRNGstate();

  for (int k = 0; k < N; ++k) {

    // rv.row(k) -= U.row(k)
    memcpy(SigmaUk, pSigmaUInv, T * T * sizeof(double));
    for (int t = 0; t < T; ++t) {
      prv[k + t * N] -= pU[k + t * N];
      omega_k[t] = pw[k + t * N];
      SigmaUk[t + t * T] += omega_k[t];
    }

    // Cholesky: SigmaUk = L L^T
    F77_CALL(dpotrf)(&uploL, &T, SigmaUk, &T, &info FCONE);
    if (info != 0) error("Cholesky failed in SigmaUk (k=%d)", k);

    // score[t] = yr[k + t*N] - w[k + t*N] * rv[k + t*N]
    for (int t = 0; t < T; ++t) {
      score[t] = pyr[k + t * N] - pw[k + t * N] * prv[k + t * N];
    }

    // Solve SigmaUk * muUk = score 
    F77_CALL(dpotrs)(&uploL, &T, &inc, SigmaUk, &T, score, &T, &info FCONE);
    if (info != 0) error("Solve linear system failed (k=%d)", k);

    // z ~ N(0, I)
    for (int t = 0; t < T; ++t) {
      z[t] = norm_rand();
    }

    // L^T y = z (i.e., z := L^{-T} z)
    F77_CALL(dtrsv)(&uploL, &transT, "N", &T, SigmaUk, &T, z, &inc FCONE FCONE FCONE);

    // rv.row(k) += U.row(k)
    for (int t = 0; t < T; ++t) {
      // U.row(k) = muUk + z
      pU[k + t * N] = score[t] + z[t];
      prv[k + t * N] += pU[k + t * N];
    }
  }

  PutRNGstate();

  SEXP result = PROTECT(allocVector(VECSXP, 2));
  SET_VECTOR_ELT(result, 0, U);
  SET_VECTOR_ELT(result, 1, rv);

  UNPROTECT(8);  // w, yr, U, SigmaUInv, rv, ntraits, n, result
  return result;
}

// Complementary functions

SEXP fllp_qmc(SEXP a_, SEXP b_, SEXP mu_, SEXP Sigma_, SEXP n_, SEXP d_, SEXP U_) {
  
  PROTECT(a_ = coerceVector(a_, REALSXP));            double *a = REAL(a_);
  PROTECT(b_ = coerceVector(b_, REALSXP));            double *b = REAL(b_);
  PROTECT(mu_ = coerceVector(mu_, REALSXP));          double *mu = REAL(mu_);
  PROTECT(Sigma_ = coerceVector(Sigma_, REALSXP));    double *Sigma = REAL(Sigma_);
  PROTECT(U_ = coerceVector(U_, REALSXP));            double *U = REAL(U_);
  PROTECT(n_ = coerceVector(n_, INTSXP));             int n = INTEGER(n_)[0];      
  PROTECT(d_ = coerceVector(d_, INTSXP));             int d = INTEGER(d_)[0];          

  int n_samples = INTEGER(getAttrib(U_, R_DimSymbol))[0]; // rows of U

  double *L = (double *) R_alloc(d * d, sizeof(double));
  memcpy(L, Sigma, d * d * sizeof(double));

  char uplo = 'L';
  int T = d;
  int info;

  F77_CALL(dpotrf)(&uplo, &T, L, &T, &info FCONE);
  if (info != 0) {
    UNPROTECT(7);
    error("Cholesky failed in Sigma (info=%d)", info);
  }

  double *z = (double *) R_alloc(d, sizeof(double));

  double loglik = 0.0;

  int i, s, j, k;

  for (i = 0; i < n; i++) {
    double acc_prob = 0.0;

    for (s = 0; s < n_samples; s++) {
      int valid = 1;
      double logp = 0.0;

      for (k = 0; k < d; k++) z[k] = 0.0;

      for (j = 0; j < d; j++) {
        double mean_j = 0.0;

        for (k = 0; k < j; k++) {
          mean_j += L[j + k * d] * z[k];
        }

        double sd_j = L[j + j * d];

        double aj = (a[i + j * n] - mu[i + j * n] - mean_j) / sd_j;
        double bj = (b[i + j * n] - mu[i + j * n] - mean_j) / sd_j;

        double cdf_a = norm_cdf(aj);
        double cdf_b = norm_cdf(bj);
        double diff = cdf_b - cdf_a;

        if (diff <= 1e-12) {
          valid = 0;
          break;
        }

        double u_val = U[s + j * n_samples];
        double u_trunc = cdf_a + u_val * diff;
        z[j] = norm_inv_cdf(u_trunc);
        logp += log(diff);
      }

      if (valid) {
        acc_prob += exp(logp);
      }
    }

    loglik += log(acc_prob / n_samples);
  }

  SEXP result = PROTECT(allocVector(REALSXP, 1));
  REAL(result)[0] = loglik;

  UNPROTECT(8);
  return result;
}

SEXP fllp_PLN_mtme(SEXP y, SEXP logy_factorial, SEXP rv, SEXP varE,
              SEXP gh_nodes, SEXP gh_weights, SEXP n_, SEXP Q_) {

  PROTECT(n_  = coerceVector(n_, INTSXP));  int n = INTEGER(n_)[0];
  PROTECT(Q_  = coerceVector(Q_, INTSXP));  int Q = INTEGER(Q_)[0];

  y             = PROTECT(coerceVector(y, REALSXP));               double* py            = REAL(y);
  logy_factorial     = PROTECT(coerceVector(logy_factorial, REALSXP));  double* plogfact      = REAL(logy_factorial);
  rv            = PROTECT(coerceVector(rv, REALSXP));              double* prv           = REAL(rv);
  gh_nodes         = PROTECT(coerceVector(gh_nodes, REALSXP));        double* pnodes        = REAL(gh_nodes);
  gh_weights       = PROTECT(coerceVector(gh_weights, REALSXP));      double* pweights      = REAL(gh_weights);

  double sigma2 = REAL(varE)[0];
  double sqrt2sigma = sqrt(2.0 * sigma2);
  double log_sqrt_pi = 0.5 * log(M_PI);
  double sumlogLik = 0.0;

  for (int i = 0; i < n; ++i) {
    double yi = py[i];
    double rvi = prv[i];

    double sum_terms = 0.0;
    for (int j = 0; j < Q; ++j) {
      double log_lambda = rvi + sqrt2sigma * pnodes[j];
      double lambda = exp(log_lambda);
      double term = exp(yi * log_lambda - lambda);
      sum_terms += pweights[j] * term;
    }

    double logLik_i = log(sum_terms) - plogfact[i] - log_sqrt_pi;
    sumlogLik += logLik_i;
  }

  UNPROTECT(7);
  return ScalarReal(sumlogLik);
}

SEXP fllp_PLN_multi(SEXP y_, SEXP logy_factorial_, SEXP rv_, 
                    SEXP gh_nodes_mult_, SEXP gh_weights_mult_,
                    SEXP n_, SEXP Q_, SEXP d_) {
 
  PROTECT(n_  = coerceVector(n_, INTSXP));  int n = INTEGER(n_)[0]; 
  PROTECT(Q_  = coerceVector(Q_, INTSXP));  int Q = INTEGER(Q_)[0];
  PROTECT(d_  = coerceVector(d_, INTSXP));  int d = INTEGER(d_)[0]; 

  SEXP y = PROTECT(coerceVector(y_, REALSXP));                                double* py = REAL(y);             
  SEXP logy_factorial = PROTECT(coerceVector(logy_factorial_, REALSXP));      double* plogfact = REAL(logy_factorial);
  SEXP rv = PROTECT(coerceVector(rv_, REALSXP));                              double* prv = REAL(rv);
  SEXP gh_nodes_mult = PROTECT(coerceVector(gh_nodes_mult_, REALSXP));        double* pnodes = REAL(gh_nodes_mult);   
  SEXP gh_weights_mult = PROTECT(coerceVector(gh_weights_mult_, REALSXP));    double* pweights = REAL(gh_weights_mult);

  // double log_sqrt_pi_d = 0.5 * d * log(M_PI);
  double sumlogLik = 0.0;

  for (int i = 0; i < n; ++i) {
    double sum_terms = 0.0;

    for (int q = 0; q < Q; ++q) {
      double dot_term = 0.0;
      for (int j = 0; j < d; ++j) {
        int idx = i + n * j;
        int qidx = q + Q * j;

        double log_lambda = prv[idx] + pnodes[qidx];
        double lambda = exp(log_lambda);

        dot_term += py[idx] * log_lambda - lambda;
      }
      sum_terms += pweights[q] * exp(dot_term);
    }
    sumlogLik += log(sum_terms) - plogfact[i];
  }

  UNPROTECT(8);
  return ScalarReal(sumlogLik);
}

SEXP fXb_multi_prod(SEXP X_, SEXP B_, SEXP rv_, SEXP alpha_) {
  double *X = REAL(X_);
  double *B = REAL(B_);
  double *rv = REAL(rv_);
  double alpha = REAL(alpha_)[0];
  double beta = 1.0;

  // Get dimensions of X
  SEXP dimX = getAttrib(X_, R_DimSymbol);
  int n = INTEGER(dimX)[0];   // number of rows
  int p = INTEGER(dimX)[1];   // number of columns

  // Get dimensions of B
  SEXP dimB = getAttrib(B_, R_DimSymbol);
  int pb = INTEGER(dimB)[0];  // number of rows
  int T = INTEGER(dimB)[1];   // number of columns
  if(pb != p) error("fXb_multi_prod: number of columns of X must equal number of rows of B");

  // Check dimensions of rv
  SEXP dimRV = getAttrib(rv_, R_DimSymbol);
  if(INTEGER(dimRV)[0] != n || INTEGER(dimRV)[1] != T)
    error("fXb_multi_prod: rv dimensions must match X * B");

  char trans = 'N';

  // Call BLAS: rv = alpha * X %*% B + beta * rv
  F77_CALL(dgemm)(
    &trans, &trans,
    &n, &T, &p,
    &alpha,
    X, &n,
    B, &p,
    &beta,
    rv, &n
    FCONE FCONE
  );

  // Return the updated rv
  return rv_;
}

SEXP prod(SEXP A_, SEXP B_) {
  // Get dimensions of A
  SEXP dimA = getAttrib(A_, R_DimSymbol);
  int m = INTEGER(dimA)[0];
  int k = INTEGER(dimA)[1];

  // Get dimensions of B
  SEXP dimB = getAttrib(B_, R_DimSymbol);
  int kb = INTEGER(dimB)[0];
  int n = INTEGER(dimB)[1];

  if (k != kb) error("Incompatible dimensions: A is %d x %d, B is %d x %d", m, k, kb, n);

  // Result object
  SEXP res = PROTECT(allocMatrix(REALSXP, m, n));
  double *C = REAL(res);

  double *A = REAL(A_);
  double *B = REAL(B_);

  char transN = 'N';
  double one = 1.0, zero = 0.0;

  // C = A %*% B
  F77_CALL(dgemm)(
    &transN, &transN,
    &m, &n, &k,
    &one,
    A, &m,
    B, &k,
    &zero,
    C, &m
    FCONE FCONE
  );

  UNPROTECT(1);
  return res;
}

SEXP kron(SEXP A_, SEXP B_) {
  // Get dimensions of A
  SEXP dimA = getAttrib(A_, R_DimSymbol);
  int m = INTEGER(dimA)[0];
  int n = INTEGER(dimA)[1];

  // Get dimensions of B
  SEXP dimB = getAttrib(B_, R_DimSymbol);
  int p = INTEGER(dimB)[0];
  int q = INTEGER(dimB)[1];

  // Allocate result matrix
  SEXP res = PROTECT(allocMatrix(REALSXP, m * p, n * q));
  double *C = REAL(res);

  const double *A = REAL(A_);
  const double *B = REAL(B_);

  int incx = 1;
  int incy = 1;

  // Loop over columns of A
  for (int j = 0; j < n; j++) {
    for (int i = 0; i < m; i++) {
      double aij = A[i + j * m];

      // Pointer to block in C (size p x q)
      double *Cblock = &C[(i * p) + (j * q) * (m * p)];

      // Copy each column of B into corresponding column of Cblock
      for (int jb = 0; jb < q; jb++) {
        const double *Bcol = &B[jb * p];
        double *Ccol = &Cblock[jb * (m * p)];

        // Copy column of B
        F77_CALL(dcopy)(&p, Bcol, &incx, Ccol, &incy);

        // Scale by aij
        F77_CALL(dscal)(&p, &aij, Ccol, &incx);
      }
    }
  }

  UNPROTECT(1);
  return res;
}

#include <R.h>
#include <Rinternals.h>

SEXP kron_vec(SEXP u_, SEXP v_) {
  int m = LENGTH(u_);
  int n = LENGTH(v_);

  const double *u = REAL(u_);
  const double *v = REAL(v_);

  SEXP res = PROTECT(allocVector(REALSXP, m * n));
  double *C = REAL(res);

  int inc = 1;

  for (int i = 0; i < m; i++) {
    double *Ci = &C[i * n];   // block for u[i]
    // copy v into block
    F77_CALL(dcopy)(&n, v, &inc, Ci, &inc);
    // scale block by u[i]
    F77_CALL(dscal)(&n, &u[i], Ci, &inc);
  }

  UNPROTECT(1);
  return res;
}