#include <stdio.h>
#include <stdbool.h>
#include <sys/param.h>
#include <math.h>
#include <memory.h>
#ifdef _OPENMP
    #include <omp.h>
#endif
#include "gsl/gsl_randist.h"
#include "gsl/gsl_cdf.h"
#include "vector.h"
#include <R.h>
#include <Rinternals.h>

#define RANGE_A 0.75
#define RANGE_B 1.5
#define RESULTS_SIZE 10
#define RESULTS_NCOL 13

double g(unsigned int s1, unsigned int r1, unsigned int n1, unsigned int s, unsigned int m, unsigned int r, unsigned int n, float p) {
    unsigned int m2 = m - n1;
    unsigned int n2 = n - n1;
    double res;
    res = gsl_cdf_binomial_P(s1, p, n1);
    for (unsigned int x = s1+1; x <= MIN(r1, s); x++) {
        res += gsl_ran_binomial_pdf(x, p, n1) * gsl_cdf_binomial_P(s-x, p, m2);
    }
    for (unsigned int x = r1+1; x <= MIN(r, n1); x++) {
        res += gsl_ran_binomial_pdf(x, p, n1) * gsl_cdf_binomial_P(r-x, p, n2);
    }
    return res;
}

SEXP g_facade(SEXP s1, SEXP r1, SEXP n1, SEXP s, SEXP m, SEXP r, SEXP n, SEXP p) {
    SEXP result = PROTECT(allocVector(REALSXP, 1));
    #pragma clang diagnostic push
    #pragma clang diagnostic ignored "-Wsign-conversion"
    REAL(result)[0] = g(asInteger(s1), asInteger(r1), asInteger(n1),
                        asInteger(s), asInteger(m),
                        asInteger(r), asInteger(n),
                        (float) asReal(p));
    #pragma clang diagnostic pop
    UNPROTECT(1);
    return result;
}

double g1(unsigned int s1, unsigned int n1, float p) {
    return gsl_cdf_binomial_P(s1, p, n1);
}

double g2(unsigned int s1, unsigned int r1, unsigned int n1, unsigned int s, unsigned int m, float p) {
    unsigned int m2 = m - n1;
    double res = 0;
    for (unsigned int x = s1+1; x <= MIN(r1, s); x++) {
        res += gsl_ran_binomial_pdf(x, p, n1) * gsl_cdf_binomial_P(s-x, p, m2);
    }
    return res;
}

double g3(unsigned int r1, unsigned int n1, unsigned int r, unsigned int n, float p) {
    unsigned int n2 = n - n1;
    double res;
    res = 0;
    for (unsigned int x = r1+1; x <= MIN(r, n1); x++) {
        res += gsl_ran_binomial_pdf(x, p, n1) * gsl_cdf_binomial_P(r-x, p, n2);
    }
    return res;
}

double en(unsigned int s1, unsigned int r1, unsigned int n1, unsigned int m2, unsigned int n2, float p) {
    return n1 +
            (gsl_cdf_binomial_P(r1, p , n1) - gsl_cdf_binomial_P(s1, p, n1)) * m2 +
            (1 - gsl_cdf_binomial_P(r1, p, n1)) * n2;
}

SEXP en_facade(SEXP s1, SEXP r1, SEXP n1, SEXP m2, SEXP n2, SEXP p) {
    SEXP result = PROTECT(allocVector(REALSXP, 1));
    #pragma clang diagnostic push
    #pragma clang diagnostic ignored "-Wsign-conversion"
    REAL(result)[0] = en(asInteger(s1), asInteger(r1), asInteger(n1),
                         asInteger(m2), asInteger(n2), (float) asReal(p));
    #pragma clang diagnostic pop
    UNPROTECT(1);
    return result;
}

typedef struct {
    unsigned int a, b;
} range;

range asmn_range(unsigned int n) {
    range r = {(unsigned int) floor(RANGE_A * n),
               (unsigned int) ceil(RANGE_B * n)};
    return r;
}

typedef struct {
    unsigned int s1, r1, n1, s, m, r, n;
    double a, b1, b2;
    double en0, en1, en2;
} params;

// does p1 better b2?
bool optimal_comparison_1(params * p1, params * p2) {
    return p1->en0 < p2->en0;
}

bool optimal_comparison_2(params * p1, params * p2) {
    return MAX(MAX(p1->en0, p1->en1), p1->en2) < MAX(MAX(p2->en0, p2->en1), p2->en2);
}

bool optimal_comparison_3(params * p1, params * p2) {
    int p1_nm = MAX(p1->n, p1->m);
    int p2_nm = MAX(p2->n, p2->m);
    return (p1_nm < p2_nm) || (p1_nm == p2_nm && p1->en0 < p2->en0);
}

bool optimal_comparison_4(params * p1, params * p2) {
    int p1_nm = MAX(p1->n, p1->m);
    int p2_nm = MAX(p2->n, p2->m);
    return (p1_nm < p2_nm) ||
        (p1_nm == p2_nm &&
            MAX(MAX(p1->en0, p1->en1), p1->en2) < MAX(MAX(p2->en0, p2->en1), p2->en2));
}

void select_optimal_results(Vector * all_params, params * opt_params[], bool (* comparison_func)(params *, params *)) {
    for (unsigned int i = 0; i < all_params->size; i++) {
        params * p = vector_get(all_params, i);
        for (unsigned int j = 0; j < RESULTS_SIZE; j++) {
            params * c = opt_params[j];
            if (c == NULL || (* comparison_func)(p, c)) {
                params * tmp = opt_params[j];
                opt_params[j] = p;
                p = tmp;
            }
        }
    }
}

void asimon(Vector * all_params, unsigned int n1, unsigned int n2, float a_max, float b1_max, float b2_max, float p0, float p1, float p2) {
    range m_range = asmn_range(n1);
    range n_range = asmn_range(n2);

    #ifdef _OPENMP
        omp_set_nested(1);
    #endif

    #pragma omp parallel for
    for (unsigned int m = m_range.a; m <= m_range.b; m++) {
        #pragma omp parallel for
        for (unsigned int n = n_range.a; n <= n_range.b; n++) {
            for (unsigned int n1 = 1; n1 <= MIN(m,n)-1; n1++) {
                for (unsigned int s1 = 0; s1 <= n1; s1++) {
                    if (s1+1 > n1) continue;
                    double b1_1 = g1(s1, n1, p1);
                    if (b1_1 > b1_max) continue;
                    double b2_1 = g1(s1, n1, p2);
                    if (b2_1 > b2_max) continue;
                    for (unsigned int r1 = s1+1; r1 <= n1; r1++) {
                        int r_bot = r1 + 1;
                        int r_up = n;
                        double a, b1, b2;

                        unsigned int r, s;
                        for (s = s1 + 1; s <= m; s++) {
                            double b1_2 = b1_1 + g2(s1, r1, n1, s, m, p1);
                            if (b1_2 > b1_max) continue;
                            double b2_2 = b2_1 + g2(s1, r1, n1, s, m, p2);
                            if (b2_2 > b2_max) continue;
                            int r_a = r_bot;
                            int r_b = r_up;
                            double b1_3, b2_3;
                            bool checked = false;
                            while (true) {
                                r = ceil((r_a + r_b) / 2.0f);
                                if (checked && r == r_up) break;
                                checked = true;
                                b1_3 = b1_2 + g3(r1, n1, r, n, p1);
                                if (b1_3 > b1_max) {
                                    r_up = r_b = r;
                                    continue;
                                }
                                b2_3 = b2_2 + g3(r1, n1, r, n, p2);
                                if (b2_3 > b2_max) {
                                    r_up = r_b = r;
                                    continue;
                                }
                                r_a = r;
                            }
                            b1 = b1_2 + g3(r1, n1, r_up - 1, n, p1);
                            b2 = b2_2 + g3(r1, n1, r_up - 1, n, p2);
                            if (b1 > b1_max || b2 > b2_max) continue;
                            if (r_up - 1 == r1) break;
                            a = 1 - g(s1, r1, n1, s, m, r_up - 1, n, p0);
                            if (a <= a_max) {
                                r_bot = r_up;
                                double en0 = en(s1, r1, n1, m - n1, n - n1, p0);
                                double en1 = en(s1, r1, n1, m - n1, n - n1, p1);
                                double en2 = en(s1, r1, n1, m - n1, n - n1, p2);
                                params params = {s1, r1, n1, s, m, r_up - 1, n,
                                                    a, b1, b2,
                                                    en0, en1, en2};
                                #pragma omp critical
                                vector_push_back(all_params, &params);
                            }
                        }
                    }
                }
            }
        }
    }
}

SEXP fill_results(params ** results_params) {
    SEXP results = PROTECT(allocMatrix(REALSXP, RESULTS_SIZE, RESULTS_NCOL));
    double * rresults = REAL(results);
    for (int i = 0; i < RESULTS_SIZE; i++) {
        params * p = results_params[i];
        rresults[i + RESULTS_SIZE * 0] = p->s1;
        rresults[i + RESULTS_SIZE * 1] = p->r1;
        rresults[i + RESULTS_SIZE * 2] = p->n1;
        rresults[i + RESULTS_SIZE * 3] = p->s;
        rresults[i + RESULTS_SIZE * 4] = p->m;
        rresults[i + RESULTS_SIZE * 5] = p->r;
        rresults[i + RESULTS_SIZE * 6] = p->n;
        rresults[i + RESULTS_SIZE * 7] = p->a;
        rresults[i + RESULTS_SIZE * 8] = p->b1;
        rresults[i + RESULTS_SIZE * 9] = p->b2;
        rresults[i + RESULTS_SIZE * 10] = p->en0;
        rresults[i + RESULTS_SIZE * 11] = p->en1;
        rresults[i + RESULTS_SIZE * 12] = p->en2;
    }
    return results;
}

SEXP asimon_facade(SEXP n1, SEXP n2, SEXP a_max, SEXP b1_max, SEXP b2_max, SEXP p0, SEXP p1, SEXP p2) {
    Vector all_params;
    vector_setup(&all_params, 1000, sizeof(params));

    #pragma clang diagnostic push
    #pragma clang diagnostic ignored "-Wsign-conversion"
    asimon(&all_params,
           asInteger(n1), asInteger(n2),
           (float) asReal(a_max), (float) asReal(b1_max), (float) asReal(b2_max),
           (float) asReal(p0), (float) asReal(p1), (float) asReal(p2));
    #pragma clang diagnostic pop

    SEXP result = PROTECT(allocVector(VECSXP, 4));

    params ** opt_params_1 = calloc(RESULTS_SIZE, sizeof(** opt_params_1));
    select_optimal_results(&all_params, opt_params_1, optimal_comparison_1);
    SET_VECTOR_ELT(result, 0, fill_results(opt_params_1));

    params ** opt_params_2 = calloc(RESULTS_SIZE, sizeof(** opt_params_2));
    select_optimal_results(&all_params, opt_params_2, optimal_comparison_2);
    SET_VECTOR_ELT(result, 1, fill_results(opt_params_2));

    params ** opt_params_3 = calloc(RESULTS_SIZE, sizeof(** opt_params_3));
    select_optimal_results(&all_params, opt_params_3, optimal_comparison_3);
    SET_VECTOR_ELT(result, 2, fill_results(opt_params_3));

    params ** opt_params_4 = calloc(RESULTS_SIZE, sizeof(** opt_params_4));
    select_optimal_results(&all_params, opt_params_4, optimal_comparison_4);
    SET_VECTOR_ELT(result, 3, fill_results(opt_params_4));

    UNPROTECT(5);

    return result;
}
