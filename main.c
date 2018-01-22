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

#define RANGE_A 0.5
#define RANGE_B 1.5
#define RESULTS_SIZE 5

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

double en(unsigned int s1, unsigned int r1, unsigned int n1, unsigned int m2, unsigned int n2, float p) {
    return n1 +
            (gsl_cdf_binomial_P(r1, p , n1) - gsl_cdf_binomial_P(s1, p, n1)) * m2 +
            (1 - gsl_cdf_binomial_P(r1, p, n1)) * n2;
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

void select_optimal_results(Vector * all_params, params ** opt_params, bool (* comparison_func)(params *, params *)) {
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

void print_results(params ** results) {
    for (int i = 0; i < RESULTS_SIZE; i++) {
        params p = * results[i];
        printf("%d %d %d | %d %d | %d %d | %.3f %.3f %.3f | %.2f %.2f %.2f\n",
               p.s1, p.r1, p.n1,
               p.s, p.m,
               p.r, p.n,
               p.a, p.b1, p.b2,
               p.en0, p.en1, p.en2);
    }
}

int main() {
    unsigned int n1, n2;
    float a_max, b1_max, b2_max;
    float p0, p1, p2;

    n1 = 29; n2 = 30;
    a_max = 0.05; b1_max = 0.2; b2_max = 0.1;
    p0 = 0.05; p1 = 0.20; p2 = 0.25;

    range m_range = asmn_range(n1);
    range n_range = asmn_range(n2);

    Vector all_params;
    vector_setup(&all_params, 1000, sizeof(params));

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
                    for(unsigned int r1 = s1+1; r1 <= n1; r1++) {
                        unsigned int i, j, ic, jc, lines, columns;
                        bool icheck, jcheck;

                        lines = n - (r1+1);
                        columns = m - (s1+1);
                        if (lines < 1 || columns < 1) continue;

                        int res[lines][columns];
                        double a;
                        double b1[lines][columns];
                        double b2[lines][columns];
                        memset(res, 0, sizeof res);

                        unsigned int r, s;
                        for (i = 0; i < lines; i++) {
                            r = r1 + 1 + i;
                            for (j = 0; j < columns; j++) {
                                s = s1 + 1 + j;

                                icheck = false;
                                for (ic = 0; ic < i; ic++) {
                                    if (res[ic][j] == -1) {
                                        icheck = true;
                                        break;
                                    }
                                }
                                if (icheck) continue;

                                jcheck = false;
                                for (jc = 0; jc < j; jc++) {
                                    if (res[i][jc] == -1) {
                                        jcheck = true;
                                        break;
                                    }
                                }
                                if (jcheck) continue;

                                b1[i][j] = g(s1, r1, n1, s, m, r, n, p1);
                                if (b1[i][j] <= b1_max) {
                                    b2[i][j] = g(s1, r1, n1, s, m, r, n, p2);
                                    if (b2[i][j] <= b2_max) {
                                        res[i][j] = 1;
                                    } else {
                                        res[i][j] = -1;
                                    }
                                } else {
                                    res[i][j] = -1;
                                }
                            }
                        }

                        for (j = 0; j < columns; j++) {
                            s = s1 + 1 + j;
                            for (i = lines-1; i != -1; i--) {
                                if (res[i][j] == 1) {
                                    r = r1 + 1 + i;
                                    a = 1 - g(s1, r1, n1, s, m, r, n, p0);
                                    if (a <= a_max) {
                                        double en0 = en(s1, r1, n1, m - n1, n - n1, p0);
                                        double en1 = en(s1, r1, n1, m - n1, n - n1, p1);
                                        double en2 = en(s1, r1, n1, m - n1, n - n1, p2);
                                        params params = {s1, r1, n1, s, m, r, n,
                                                         a, b1[i][j], b2[i][j],
                                                         en0, en1, en2};
                                        #pragma omp critical
                                        vector_push_back(&all_params, &params);
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    params * opt_params_1[RESULTS_SIZE];
    select_optimal_results(&all_params, opt_params_1, optimal_comparison_1);
    print_results(opt_params_1);

    printf("\n");

    params * opt_params_2[RESULTS_SIZE];
    select_optimal_results(&all_params, opt_params_2, optimal_comparison_2);
    print_results(opt_params_2);

    printf("\n");

    params * opt_params_3[RESULTS_SIZE];
    select_optimal_results(&all_params, opt_params_3, optimal_comparison_3);
    print_results(opt_params_3);

    printf("\n");

    params * opt_params_4[RESULTS_SIZE];
    select_optimal_results(&all_params, opt_params_4, optimal_comparison_4);
    print_results(opt_params_4);

    return 0;
}
