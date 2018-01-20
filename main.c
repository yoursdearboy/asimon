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

bool check_b1(unsigned int s1, unsigned int r1, unsigned int n1, unsigned int s, unsigned int m, unsigned int r, unsigned int n, float p1, float b1) {
    if (g(s1, r1, n1, s, m, r, n, p1) <= b1) {
        return true;
    } else {
        return false;
    }
}

bool check_b2(unsigned int s1, unsigned int r1, unsigned int n1, unsigned int s, unsigned int m, unsigned int r, unsigned int n, float p2, float b2) {
    if (g(s1, r1, n1, s, m, r, n, p2) <= b2) {
        return true;
    } else {
        return false;
    }
}

bool check_a(unsigned int s1, unsigned int r1, unsigned int n1, unsigned int s, unsigned int m, unsigned int r, unsigned int n, float p0, float a) {
    if (g(s1, r1, n1, s, m, r, n, p0) > (1-a)) {
        return true;
    } else {
        return false;
    }
}

typedef struct {
    unsigned int a, b;
} range;

range asmn_range(unsigned int n) {
    range r = {(unsigned int) floor(0.85 * n),
               (unsigned int) ceil(1.5 * n)};
    return r;
}

typedef struct {
    unsigned int s1, r1, n1, s, m, r, n;
} params;

int main() {
    unsigned int n1, n2;
    float a, b1, b2;
    float p0, p1, p2;

    n1 = 29; n2 = 41;
    a = 0.05; b1 = 0.2; b2 = 0.1;
    p0 = 0.05; p1 = 0.20; p2 = 0.25;

    range n_range = asmn_range(n1);
    range m_range = asmn_range(n2);

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

                                if (check_b1(s1, r1, n1, s, m, r, n, p1, b1) &&
                                    check_b2(s1, r1, n1, s, m, r, n, p2, b2)) {
                                    res[i][j] = 1;
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
                                    if (check_a(s1, r1, n1, s, m, r, n, p0, a)) {
                                        params params = {s1, r1, n1, s, m, r, n};
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

    VECTOR_FOR_EACH(&all_params, i) {
        params p = ITERATOR_GET_AS(params, &i);
        printf("%d %d %d %d %d\n", p.s1, p.r1, p.n1, p.s, p.r);
    }

    return 0;
}
