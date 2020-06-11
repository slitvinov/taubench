#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "types.h"
#include "main.h"
#include "smooth.h"
#include "flop.h"
#include "util.h"
void
kernel_3_0(DualGrid * grid, WorkSpace * work, double dummy_4[][CO],
           double dummy_3[][AD], double dummy_1[][PR], int sm, int re,
           int id, int prn)
{
    RangeList *cl;
    double *fin_c = NULL, *fin_v = NULL;
    const double four_thirds = FOUR_THIRDS;
    double (*hop)[3] = grid->hop;
    double (*wfl)[3] = work->wfl;
    int (*hup)[2] = grid->hup;
    int viscous = 1;
    int moving_grid = grid->moving;
    int skip;
    const int npnt = grid->nallpoints;
    int nskips = grid->nallskips;
    int step;
    static double st1 = 0.0;
    static int stc = 0;
    double t1 = 0.0, t2 = 0.0;
    double junk_1 = 1.0;
    double junk_2 = 1.0;
    double junk_3 = 1.0;
    double junk_4 = 1.0;
    double junk_5 = 1.0;

    if (id == 0) {
        stc++;
        t1 = second();
    }
    fin_c = (double (*)) test_malloc(nskips * sizeof(double));
    fin_v = (double (*)) test_malloc(nskips * sizeof(double));
    for (step = 0; step < re; step++) {
        for (cl = grid->fcl; cl != NULL; cl = cl->succ) {
            const int start = cl->start, stop = cl->stop;

#include "nodep.h"
            for (skip = start; skip < stop; skip++) {
                const int p0 = hup[skip][0], p1 = hup[skip][1];
                const double anx = hop[skip][0], any = hop[skip][1], anz =
                    hop[skip][2];
                const double area = BLA_TWO(hop, skip);
                const double a0 = BLA_SEVEN(dummy_1, p0), a1 =
                    BLA_SEVEN(dummy_1, p1);
                double v0 =
                    dummy_1[p0][IVX] * anx + dummy_1[p0][IVY] * any +
                    dummy_1[p0][IVZ] * anz, v1 =
                    dummy_1[p1][IVX] * anx + dummy_1[p1][IVY] * any +
                    dummy_1[p1][IVZ] * anz;
                double a = 0.5 * (a0 + a1);
                double vm = 0.5 * (v0 + v1);

                if (moving_grid) {
                    const double wfln =
                        0.5 * ((wfl[p0][0] + wfl[p1][0]) * anx +
                               (wfl[p0][1] + wfl[p1][1]) * any +
                               (wfl[p0][2] + wfl[p1][2]) * anz);
                    vm -= wfln;
                }
                fin_c[skip] = fabs(vm) + area * a;
                if (viscous) {
                    const double rho =
                        0.5 * (dummy_1[p0][IRHO] + dummy_1[p1][IRHO]);
                    const double tem0 = BLA_THREE(dummy_1, p0), tem1 =
                        BLA_THREE(dummy_1, p1);
                    const double mue_l0 =
                        BLA_ONE(tem0, dummy_1, p0), mue_l1 =
                        BLA_ONE(tem1, dummy_1, p1), kappa_l0 =
                        BLA_FOUR(tem0, dummy_1, p0), kappa_l1 =
                        BLA_FOUR(tem1, dummy_1, p1), c_v0 =
                        BLA_SIX(dummy_1, p0), c_v1 = BLA_SIX(dummy_1, p1);
                    double mue_lam = 0.5 * (mue_l0 + mue_l1), kappa_lam =
                        0.5 * (kappa_l0 + kappa_l1), c_v =
                        0.5 * (c_v0 + c_v1);
                    double fin_v1 = mue_lam, fin_v2 = kappa_lam, fin_max;
                    const double eddyv =
                        0.5 * (dummy_3[p0][AEDV] + dummy_3[p1][AEDV]);
                    fin_v1 += eddyv;
                    fin_v2 +=
                        grid->ptl_ratio * kappa_lam / mue_lam * eddyv;
                    fin_v1 *= four_thirds;
                    fin_v2 /= c_v;
                    fin_max = fin_v1 + fin_v2;
                    fin_v[skip] = area * area * fin_max / rho;
                }
    }}} free(fin_v);
    free(fin_c);
    if (id == 0) {
        t2 = second();
        st1 += (t2 - t1);
        if (prn != 0) {
            printf("        - kernel_3_0 : %10.3f secs - %10.3f mflops\n",
                   st1, SAL * npnt / st1 * stc);
        }
    }
}
