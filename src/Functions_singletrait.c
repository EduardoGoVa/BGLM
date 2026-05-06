#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/BLAS.h> 
#include <R_ext/Lapack.h>
#include <R_ext/Rdynload.h>

// component-component update functions

SEXP fBj_DLN_G(SEXP Bv, SEXP varB, SEXP varE, SEXP rv, SEXP X, SEXP x2, 
           SEXP n, SEXP p) {

  int rows = INTEGER(n)[0];  
  int cols = INTEGER(p)[0];
  int inc = 1;

  double svarE;
  double Bk, rhs, c, diff;

  double* pBv;
  double* pvarB;
  double* prv;
  double* pX;
  double* px2;
  
  PROTECT(Bv = coerceVector(Bv, REALSXP));        pBv = REAL(Bv);
  PROTECT(varB = coerceVector(varB, REALSXP));    pvarB = REAL(varB);
  PROTECT(rv = coerceVector(rv, REALSXP));        prv = REAL(rv);
  PROTECT(X = coerceVector(X, REALSXP));          pX = REAL(X);
  PROTECT(x2 = coerceVector(x2, REALSXP));        px2 = REAL(x2);

  svarE = 1.0 / REAL(varE)[0];

  GetRNGstate();

  for (int k = 0; k < cols; ++k) {
    double* Xk = pX + (long long)k * rows;
    Bk = pBv[k];

    rhs = F77_CALL(ddot)(&rows, prv, &inc, Xk, &inc) * svarE;
    
    rhs += px2[k] * Bk * svarE;

    c = px2[k] * svarE + 1.0 / pvarB[k];

    pBv[k] = rhs / c + sqrt(1.0 / c) * norm_rand();

    diff = Bk - pBv[k];

    F77_CALL(daxpy)(&rows, &diff, Xk, &inc, prv, &inc);
  }

  PutRNGstate();

  SEXP res = PROTECT(allocVector(VECSXP, 2));
  SET_VECTOR_ELT(res, 0, Bv);
  SET_VECTOR_ELT(res, 1, rv);

  UNPROTECT(6);
  return res;
}

SEXP fUi_DLN_G(SEXP Bv, SEXP varB, SEXP varE, SEXP rv, SEXP X, SEXP n, SEXP p) {

  int rows = INTEGER(n)[0];  
  int cols = INTEGER(p)[0];
  int inc = 1;

  double svarE;
  double Bk, rhs, c, diff;

  double* pBv;
  double* pvarB;
  double* prv;
  double* pX;
  
  PROTECT(Bv = coerceVector(Bv, REALSXP));        pBv = REAL(Bv);
  PROTECT(varB = coerceVector(varB, REALSXP));    pvarB = REAL(varB);
  PROTECT(rv = coerceVector(rv, REALSXP));        prv = REAL(rv);
  PROTECT(X = coerceVector(X, REALSXP));          pX = REAL(X);

  svarE = 1.0 / REAL(varE)[0];

  GetRNGstate();

  for (int k = 0; k < cols; ++k) {
    double* Xk = pX + (long long)k * rows;
    Bk = pBv[k];

    rhs = F77_CALL(ddot)(&rows, prv, &inc, Xk, &inc) * svarE;
    
    rhs += Bk * svarE;

    c = svarE + 1.0 / pvarB[k];

    pBv[k] = rhs / c + sqrt(1.0 / c) * norm_rand();

    diff = Bk - pBv[k];

    F77_CALL(daxpy)(&rows, &diff, Xk, &inc, prv, &inc);
  }

  PutRNGstate();

  SEXP res = PROTECT(allocVector(VECSXP, 2));
  SET_VECTOR_ELT(res, 0, Bv);
  SET_VECTOR_ELT(res, 1, rv);

  UNPROTECT(5);
  return res;
}

SEXP fBj_P(SEXP w, SEXP yr, SEXP Bv, SEXP varB, SEXP rv, SEXP X, SEXP X2, 
           SEXP n, SEXP p) {
  int rows = INTEGER(n)[0];
  int cols = INTEGER(p)[0];
  int inc = 1;

  double mu_k, var_k, Bk, diff;

  double *pw, *pyr, *pBv, *pvarB, *prv, *pX, *pX2;

  PROTECT(w       = coerceVector(w,        REALSXP)); pw       = REAL(w);
  PROTECT(yr      = coerceVector(yr,       REALSXP)); pyr      = REAL(yr);
  PROTECT(Bv      = coerceVector(Bv,       REALSXP)); pBv      = REAL(Bv);
  PROTECT(varB    = coerceVector(varB,     REALSXP)); pvarB    = REAL(varB);
  PROTECT(rv      = coerceVector(rv,       REALSXP)); prv      = REAL(rv);
  PROTECT(X       = coerceVector(X,        REALSXP)); pX       = REAL(X);
  PROTECT(X2      = coerceVector(X2,       REALSXP)); pX2      = REAL(X2);

  GetRNGstate();

  for (int k = 0; k < cols; ++k) {
    double* Xk  = pX  + (long long)k * rows;
    double* X2k = pX2 + (long long)k * rows;

    Bk = -pBv[k];

    double dot_val = 0.0;
    for (int i = 0; i < rows; ++i) {
      dot_val += pw[i] * (Xk[i] * prv[i] + X2k[i] * Bk);
    }

    double dot_w   = F77_CALL(ddot)(&rows, pw,  &inc, X2k,    &inc);
    double dot_yr  = F77_CALL(ddot)(&rows, pyr, &inc, Xk, &inc);

    var_k = 1.0 / (1.0 / pvarB[k] + dot_w);
    mu_k  = var_k * (dot_yr - dot_val);

    pBv[k] = mu_k + sqrt(var_k) * norm_rand();

    diff = Bk + pBv[k];

    F77_CALL(daxpy)(&rows, &diff, Xk, &inc, prv, &inc);  // rv += Bk * Xk
  }

  PutRNGstate();

  SEXP res = PROTECT(allocVector(VECSXP, 2));
  SET_VECTOR_ELT(res, 0, Bv);
  SET_VECTOR_ELT(res, 1, rv);

  UNPROTECT(8); 

  return res;
}

SEXP fUi(SEXP w, SEXP yr, SEXP U, SEXP varE, SEXP rv, SEXP n_) {
  int n = INTEGER(n_)[0];
  double svarU = 1.0 / REAL(varE)[0];

  double *pw, *pyr, *pU, *prv;
  PROTECT(w  = coerceVector(w,  REALSXP)); pw = REAL(w);
  PROTECT(yr = coerceVector(yr, REALSXP)); pyr = REAL(yr);
  PROTECT(U  = coerceVector(U,  REALSXP)); pU = REAL(U);
  PROTECT(rv = coerceVector(rv, REALSXP)); prv = REAL(rv);

  GetRNGstate();

  for (int k = 0; k < n; ++k) {
    prv[k] -= pU[k];

    double sigma2_Uk = 1.0 / (svarU + pw[k]);
    double mu_Uk = sigma2_Uk * (pyr[k] - pw[k] * prv[k]);

    pU[k] = mu_Uk + sqrt(sigma2_Uk) * norm_rand();

    prv[k] += pU[k];
  }

  PutRNGstate();

  SEXP res = PROTECT(allocVector(VECSXP, 2));
  SET_VECTOR_ELT(res, 0, U);
  SET_VECTOR_ELT(res, 1, rv);

  UNPROTECT(5); // 4 PROTECTs + res
  return res;
}

// block update functions

SEXP fB_DLN_G(SEXP Bv, SEXP SigmaBInv, SEXP rv,
            SEXP X, SEXP tXX, SEXP tX, SEXP varE, SEXP blocks,
            SEXP nblocks, SEXP n) {

  PROTECT(Bv        = coerceVector(Bv, REALSXP));
  PROTECT(SigmaBInv = coerceVector(SigmaBInv, REALSXP));
  PROTECT(rv        = coerceVector(rv, REALSXP));
  PROTECT(X         = coerceVector(X, REALSXP));
  PROTECT(tXX       = coerceVector(tXX, REALSXP));
  PROTECT(tX        = coerceVector(tX, REALSXP));
  PROTECT(varE      = coerceVector(varE, REALSXP));
  PROTECT(nblocks   = coerceVector(nblocks, INTSXP));
  PROTECT(n         = coerceVector(n, INTSXP));

  const int N = INTEGER(n)[0];              // Total number of observations
  const int K = INTEGER(nblocks)[0];
  const int P = LENGTH(Bv);

  const double svarE = 1.0 / REAL(varE)[0];

  double *pBv = REAL(Bv);
  double *pSigmaBInv = REAL(SigmaBInv);
  double *prv = REAL(rv);
  double *pX = REAL(X);
  double *ptXX = REAL(tXX);
  double *ptX = REAL(tX);

  const int inc = 1;
  int info;
  const double one = 1.0, zero = 0.0, neg_one = -1.0;
  const char Trans = 'N', TransT = 'T', Uplo = 'L';

  SEXP result;
  PROTECT(result = allocVector(VECSXP, 2));

  GetRNGstate();

  for (int k = 0; k < K; ++k) {
    SEXP idx_block = VECTOR_ELT(blocks, k);
    const int m = LENGTH(idx_block);
    int* pidx = INTEGER(idx_block);

    double *Bv_block = (double*) R_alloc(m, sizeof(double));
    double *X_block = (double*) R_alloc(N * m, sizeof(double));
    double *tX_block = (double*) R_alloc(m * N, sizeof(double));  
    double *S = (double*) R_alloc(m * m, sizeof(double));
    double *mu = (double*) R_alloc(m, sizeof(double));
    double *z = (double*) R_alloc(m, sizeof(double));
    double *Bk = (double*) R_alloc(m, sizeof(double));           

    // Build Bv_block, X_block, tX_block and matrix S
    for (int j = 0; j < m; ++j) {
      int col_j = pidx[j] - 1;
      Bv_block[j] = pBv[col_j];

      // Copy column col_j of X for all observations
      for (int i = 0; i < N; ++i) {
        X_block[i + j * N] = pX[i + col_j * N];
        tX_block[j + i * m] = ptX[col_j + i * P];
      }

      // Build matrix S = SigmaBInv + svarE * tXX_block (m x m)
      for (int l = 0; l < m; ++l) {
        int col_l = pidx[l] - 1;
        S[j + l * m] = pSigmaBInv[col_j + col_l * P] + svarE * ptXX[col_j + col_l * P];
      }
    }

    // rv += X_block * Bv_block
    F77_CALL(dgemv)(&Trans, &N, &m, &one, X_block, &N, Bv_block, &inc, &one, prv, &inc FCONE);

    // Cholesky decomposition: S = L * Lᵀ
    F77_CALL(dpotrf)(&Uplo, &m, S, &m, &info FCONE);
    if (info != 0)
      error("Cholesky failed at block %d", k + 1);

    // mu = solve(S, svarE * tX_block %*% lr)
    F77_CALL(dgemv)(&Trans, &m, &N, &svarE, tX_block, &m, prv, &inc, &zero, mu, &inc FCONE);
    F77_CALL(dpotrs)(&Uplo, &m, &inc, S, &m, mu, &m, &info FCONE);
    if (info != 0)
      error("Solve failed at block %d", k + 1);

    // z ~ N(0, I), solve Lᵀ z = z
    for (int i = 0; i < m; ++i)
      z[i] = norm_rand();
    F77_CALL(dtrsv)(&Uplo, &TransT, "N", &m, S, &m, z, &inc FCONE FCONE FCONE);

    // Update Bv block: Bk = mu + z
    for (int i = 0; i < m; ++i) {
      Bk[i] = mu[i] + z[i];
      int col = pidx[i] - 1;
      pBv[col] = Bk[i];
    }

    // rv -= X_block * Bk
    F77_CALL(dgemv)(&Trans, &N, &m, &neg_one, X_block, &N, Bk, &inc, &one, prv, &inc FCONE);
  }

  PutRNGstate();

  SET_VECTOR_ELT(result, 0, Bv);
  SET_VECTOR_ELT(result, 1, rv);

  UNPROTECT(10);
  return result;
}

SEXP fU_DLN_G(SEXP Bv, SEXP SigmaBInv, SEXP rv,
            SEXP X, SEXP tX, SEXP varE, SEXP blocks,
            SEXP nblocks, SEXP n) {

  PROTECT(Bv        = coerceVector(Bv, REALSXP));
  PROTECT(SigmaBInv = coerceVector(SigmaBInv, REALSXP));
  PROTECT(rv        = coerceVector(rv, REALSXP));
  PROTECT(X         = coerceVector(X, REALSXP));
  PROTECT(tX        = coerceVector(tX, REALSXP));
  PROTECT(varE      = coerceVector(varE, REALSXP));
  PROTECT(nblocks   = coerceVector(nblocks, INTSXP));
  PROTECT(n         = coerceVector(n, INTSXP));

  const int N = INTEGER(n)[0];              // Total number of observations
  const int K = INTEGER(nblocks)[0];
  const int P = LENGTH(Bv);

  const double svarE = 1.0 / REAL(varE)[0];

  double *pBv = REAL(Bv);
  double *pSigmaBInv = REAL(SigmaBInv);
  double *prv = REAL(rv);
  double *pX = REAL(X);
  double *ptX = REAL(tX);

  const int inc = 1;
  int info;
  const double one = 1.0, zero = 0.0, neg_one = -1.0;
  const char Trans = 'N', TransT = 'T', Uplo = 'L';

  SEXP result;
  PROTECT(result = allocVector(VECSXP, 2));

  GetRNGstate();

  for (int k = 0; k < K; ++k) {
    SEXP idx_block = VECTOR_ELT(blocks, k);
    const int m = LENGTH(idx_block);
    int* pidx = INTEGER(idx_block);

    double *Bv_block = (double*) R_alloc(m, sizeof(double));
    double *X_block = (double*) R_alloc(N * m, sizeof(double));
    double *tX_block = (double*) R_alloc(m * N, sizeof(double));  
    double *S = (double*) R_alloc(m * m, sizeof(double));
    double *mu = (double*) R_alloc(m, sizeof(double));
    double *z = (double*) R_alloc(m, sizeof(double));
    double *Bk = (double*) R_alloc(m, sizeof(double));           

    // Build Bv_block, X_block, tX_block and matrix S
    for (int j = 0; j < m; ++j) {
      int col_j = pidx[j] - 1;
      Bv_block[j] = pBv[col_j];

      // Copy column col_j of X for all observations
      for (int i = 0; i < N; ++i) {
        X_block[i + j * N] = pX[i + col_j * N];
        tX_block[j + i * m] = ptX[col_j + i * P];
      }

      // Build matrix S = SigmaBInv + svarE
      for (int l = 0; l < m; ++l) {
        int col_l = pidx[l] - 1;
        S[j + l * m] = pSigmaBInv[col_j + col_l * P] + svarE;
      }
    }

    // rv += X_block * Bv_block
    F77_CALL(dgemv)(&Trans, &N, &m, &one, X_block, &N, Bv_block, &inc, &one, prv, &inc FCONE);

    // Cholesky decomposition: S = L * Lᵀ
    F77_CALL(dpotrf)(&Uplo, &m, S, &m, &info FCONE);
    if (info != 0)
      error("Cholesky failed at block %d", k + 1);

    // mu = solve(S, svarE * tX_block %*% lr)
    F77_CALL(dgemv)(&Trans, &m, &N, &svarE, tX_block, &m, prv, &inc, &zero, mu, &inc FCONE);
    F77_CALL(dpotrs)(&Uplo, &m, &inc, S, &m, mu, &m, &info FCONE);
    if (info != 0)
      error("Solve failed at block %d", k + 1);

    // z ~ N(0, I), solve Lᵀ z = z
    for (int i = 0; i < m; ++i)
      z[i] = norm_rand();
    F77_CALL(dtrsv)(&Uplo, &TransT, "N", &m, S, &m, z, &inc FCONE FCONE FCONE);

    // Update Bv block: Bk = mu + z
    for (int i = 0; i < m; ++i) {
      Bk[i] = mu[i] + z[i];
      int col = pidx[i] - 1;
      pBv[col] = Bk[i];
    }

    // rv -= X_block * Bk
    F77_CALL(dgemv)(&Trans, &N, &m, &neg_one, X_block, &N, Bk, &inc, &one, prv, &inc FCONE);
  }

  PutRNGstate();

  SET_VECTOR_ELT(result, 0, Bv);
  SET_VECTOR_ELT(result, 1, rv);

  UNPROTECT(9);
  return result;
}

SEXP fB_P(SEXP w, SEXP yr, SEXP Bv, SEXP SigmaBInv, 
          SEXP rv, SEXP X, SEXP tX, SEXP blocks, 
          SEXP nblocks, SEXP n) {

  PROTECT(w = coerceVector(w, REALSXP));
  PROTECT(yr = coerceVector(yr, REALSXP));
  PROTECT(Bv = coerceVector(Bv, REALSXP));
  PROTECT(SigmaBInv = coerceVector(SigmaBInv, REALSXP));
  PROTECT(rv = coerceVector(rv, REALSXP));
  PROTECT(X = coerceVector(X, REALSXP));
  PROTECT(tX = coerceVector(tX, REALSXP));
  PROTECT(nblocks = coerceVector(nblocks, INTSXP));
  PROTECT(n = coerceVector(n, INTSXP));

  int N = INTEGER(n)[0];
  int K = INTEGER(nblocks)[0];
  int P = LENGTH(Bv);

  double *pw = REAL(w);
  double *pyr = REAL(yr);
  double *pBv = REAL(Bv);
  double *pSigmaBInv = REAL(SigmaBInv);
  double *prv = REAL(rv);
  double *pX = REAL(X);
  double *ptX = REAL(tX);

  int inc = 1;
  double one = 1.0, zero = 0.0, neg_one = -1.0;
  double alpha = 1.0, beta = 1.0;
  char trans = 'N', uplo = 'L';
  int info;

  SEXP result;
  PROTECT(result = allocVector(VECSXP, 2));

  GetRNGstate();

  for (int k = 0; k < K; ++k) {
    SEXP idx = VECTOR_ELT(blocks, k);
    int m = LENGTH(idx);
    int* pidx = INTEGER(idx);

    double* Bv_block = (double*) R_alloc(m, sizeof(double));
    double* X_block = (double*) R_alloc(N * m, sizeof(double));
    double* tX_block = (double*) R_alloc(m * N, sizeof(double));
    double* Xw = (double*) R_alloc(N * m, sizeof(double));
    double* DtXwt = (double*) R_alloc(m * N, sizeof(double));
    double* S = (double*) R_alloc(m * m, sizeof(double));
    double* temp1 = (double*) R_alloc(m, sizeof(double));
    double* temp2 = (double*) R_alloc(m, sizeof(double));
    double* mu = (double*) R_alloc(m, sizeof(double));
    double* z = (double*) R_alloc(m, sizeof(double));
    double* Bk = (double*) R_alloc(m, sizeof(double));
    double* SigmaBInv_block = (double*) R_alloc(m * m, sizeof(double));

    // Build Bv_block, X_block, tX_block and S
    for (int j = 0; j < m; ++j) {
      int col_j = pidx[j] - 1;
      Bv_block[j] = pBv[col_j];

      for (int i = 0; i < N; ++i) {
        X_block[i + j * N] = pX[i + col_j * N];
        tX_block[j + i * m] = ptX[col_j + i * P];
        double sqw = sqrt(pw[i]);
        Xw[i + j * N] = sqw * pX[i + col_j * N];
        DtXwt[j + i * m] = pw[i] * pX[i + col_j * N];
      }

      for (int l = 0; l < m; ++l) {
        int col_l = pidx[l] - 1;
        SigmaBInv_block[j + l * m] = pSigmaBInv[col_j + col_l * P];
      }
    }

    // rv -= X_block * Bv_block
    F77_CALL(dgemv)(&trans, &N, &m, &neg_one, X_block, &N, Bv_block, &inc, &one, prv, &inc FCONE);

    // S = SigmaBInv_block + Xw' * Xw using dgemm
    memcpy(S, SigmaBInv_block, m * m * sizeof(double));
    F77_CALL(dgemm)("T", "N", &m, &m, &N, &alpha, Xw, &N, Xw, &N, &beta, S, &m FCONE FCONE);

    // Cholesky decomposition
    F77_CALL(dpotrf)(&uplo, &m, S, &m, &info FCONE);
    if (info != 0)
      error("Cholesky failed at block %d", k + 1);

    // temp1 = tX_block * yr
    F77_CALL(dgemv)(&trans, &m, &N, &one, tX_block, &m, pyr, &inc, &zero, temp1, &inc FCONE);

    // temp2 = DtXwt * rv
    F77_CALL(dgemv)(&trans, &m, &N, &one, DtXwt, &m, prv, &inc, &zero, temp2, &inc FCONE);

    for (int i = 0; i < m; ++i)
      mu[i] = temp1[i] - temp2[i];

    // Solve S * mu = mu (using Cholesky)
    F77_CALL(dpotrs)(&uplo, &m, &inc, S, &m, mu, &m, &info FCONE);
    if (info != 0)
      error("Solve failed at block %d", k + 1);

    // Generate z ~ N(0, I)
    for (int i = 0; i < m; ++i)
      z[i] = norm_rand();

    // Solve L^T z = z (where L is Cholesky factor)
    F77_CALL(dtrsv)(&uplo, "T", "N", &m, S, &m, z, &inc FCONE FCONE FCONE);

    // Update Bk = mu + z and assign back to Bv
    for (int i = 0; i < m; ++i) {
      Bk[i] = mu[i] + z[i];
      pBv[pidx[i] - 1] = Bk[i];
    }

    // rv += X_block * Bk
    F77_CALL(dgemv)(&trans, &N, &m, &one, X_block, &N, Bk, &inc, &one, prv, &inc FCONE);
  }

  PutRNGstate();

  SET_VECTOR_ELT(result, 0, Bv);
  SET_VECTOR_ELT(result, 1, rv);

  UNPROTECT(10);
  return result;
}

// Complementary functions

SEXP fllp_PLN(SEXP y, SEXP logy_factorial, SEXP rv, SEXP varE,
              SEXP gh_nodes, SEXP gh_weights, SEXP n_, SEXP Q_) {
  int n = INTEGER(n_)[0];
  int Q = INTEGER(Q_)[0];

  SEXP vecy             = PROTECT(coerceVector(y, REALSXP));
  SEXP veclogy_fact     = PROTECT(coerceVector(logy_factorial, REALSXP));
  SEXP vecrv            = PROTECT(coerceVector(rv, REALSXP));
  SEXP vecnodes         = PROTECT(coerceVector(gh_nodes, REALSXP));
  SEXP vecweights       = PROTECT(coerceVector(gh_weights, REALSXP));

  double* py            = REAL(vecy);
  double* plogfact      = REAL(veclogy_fact);
  double* prv           = REAL(vecrv);
  double* pnodes        = REAL(vecnodes);
  double* pweights      = REAL(vecweights);

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

  UNPROTECT(5);
  return ScalarReal(sumlogLik);
}

SEXP fXb_prod(SEXP X_, SEXP B_, SEXP mu_, SEXP n_, SEXP p_) {
  double *X = REAL(X_);
  double *B = REAL(B_);
  double *mu = REAL(mu_);
  
  int n = INTEGER(n_)[0];
  int p = INTEGER(p_)[0];
  int inc = 1;
  char trans = 'N';
  double alpha = 1.0;
  double beta = 1.0;

  F77_CALL(dgemv)(&trans, &n, &p, &alpha, X, &n, B, &inc, &beta, mu, &inc FCONE);

  return R_NilValue;
}
