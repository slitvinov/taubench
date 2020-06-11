#include <math.h>
#include "malloc.h"
#include "util.h"
#include "types.h"
#include "main.h"
#include "flux.h"
#include "flop.h"
static void kernel_2_1(DualGrid * grid, WorkSpace * work,
                       double dummy_4[][CO]);
static void kernel_2_2(DualGrid * grid, double dummy_1[][PR],
                       double dummy_3[][AD], double dummy_2[][GR][3],
                       double dummy_4[][CO]);
static void kernel_2_3(DualGrid * grid, double dummy_1[][PR],
                       double dummy_3[][AD], double dummy_2[][GR][3],
                       double dummy_4[][CO]);
static void kernel_2_4(DualGrid * grid, double dummy_1[][PR],
                       double dummy_3[][AD], double dummy_2[][GR][3],
                       double dummy_4[][CO]);
void
kernel_2_0(DualGrid * grid, WorkSpace * work, double dummy_1[][PR],
           double dummy_3[][AD], double dummy_2[][GR][3],
           double dummy_4[][CO], int id, int prn)
{
    const int npnt = grid->nallpoints;
    double t1 = 0.0, t2 = 0.0;
    static double st1 = 0.0, st2 = 0.0, st3 = 0.0, st4 = 0.0;
    static int stc = 0;

    if (id == 0) {
        stc++;
        t1 = second();
    }
    kernel_2_1(grid, work, dummy_4);
    if (id == 0) {
        t2 = second();
        st1 += (t2 - t1);
        if (prn != 0) {
            printf("        - kernel_2_1 : %10.3f secs - %10.3f mflops\n",
                   st1, ZON * npnt / st1 * stc);
        }
    }
    if (id == 0) {
        t1 = second();
    }
    kernel_2_2(grid, dummy_1, dummy_3, dummy_2, dummy_4);
    if (id == 0) {
        t2 = second();
        st2 += (t2 - t1);
        if (prn != 0) {
            printf("        - kernel_2_2 : %10.3f secs - %10.3f mflops\n",
                   st2, ZTW * npnt / st2 * stc);
        }
    }
    if (id == 0) {
        t1 = second();
    }
    kernel_2_3(grid, dummy_1, dummy_3, dummy_2, dummy_4);
    if (id == 0) {
        t2 = second();
        st3 += (t2 - t1);
        if (prn != 0) {
            printf("        - kernel_2_3 : %10.3f secs - %10.3f mflops\n",
                   st3, ZTH * npnt / st3 * stc);
        }
    }
    if (id == 0) {
        t1 = second();
    }
    kernel_2_4(grid, dummy_1, dummy_3, dummy_2, dummy_4);
    if (id == 0) {
        t2 = second();
        st4 += (t2 - t1);
        if (prn != 0) {
            printf("        - kernel_2_4 : %10.3f secs - %10.3f mflops\n",
                   st4, ZFO * npnt / st4 * stc);
        }
    }
}

static void
kernel_2_1(DualGrid * grid, WorkSpace * work, double dummy_4[][CO])
{
    double (*hop)[3] = grid->hop;
    double (*wfl)[3] = work->wfl;
    const int moving_grid = grid->moving;
    int (*hup)[2] = grid->hup;
    double (*ping)[ST], (*pong)[ST];
    RangeList *cl;
    int nskips = grid->nallskips;

    ping = (double (*)[ST]) test_malloc(nskips * sizeof(double[ST]));
    pong = (double (*)[ST]) test_malloc(nskips * sizeof(double[ST]));
    for (cl = grid->fcl; cl != NULL; cl = cl->succ) {
        const int start = cl->start;
        const int stop = cl->stop;
        int skip, j;

        for (skip = start; skip < stop; skip++) {
#define EXPAND_COUNT ST
#include "expand.h"
            for (j = 0; j < ST; j++) {
                ping[skip][j] = 1.0;
                pong[skip][j] = 1.0;
            }
        }
    }
    for (cl = grid->fcl; cl != NULL; cl = cl->succ) {
        const int start = cl->start;
        const int stop = cl->stop;
        int skip;

#include "nodep.h"
        for (skip = start; skip < stop; skip++) {
            int k;
            const int pl = hup[skip][0], pr = hup[skip][1];
            const double anx = hop[skip][0], any = hop[skip][1], anz =
                hop[skip][2];
            const double area = BLA_TWO(hop, skip);
            const double rarea = 1 / area;
            const double nx = anx * rarea, ny = any * rarea, nz =
                anz * rarea;
            double rho_l, vx_l, vy_l, vz_l, p_l, a_l, m_l, rho_r, vx_r,
                vy_r, vz_r, p_r, a_r, m_r, p_12 = 0.0, tmp, ssw =
                0.0, ssw_a, rhovn_l = 0.0, rhovn_r =
                0.0, p_12_l, p_12_r, vn_l, vn_r, gs =
                0.0, vn_gs_l, vn_gs_r, a_m, vn_pl, vn_mi, p_pl, p_mi,
                rhovn, rhovvn = 0.0, s;
            double rhovn_la, rhovn_ra;
            double zap[CO], var_r[ZAP], var_l[ZAP];
            int atmp, btmp;

            rho_l = ping[skip][IRHO];
            vx_l = ping[skip][IVX];
            vy_l = ping[skip][IVY];
            vz_l = ping[skip][IVZ];
            p_l = ping[skip][IP];
            a_l = ping[skip][IA];
            vn_l = vx_l * nx + vy_l * ny + vz_l * nz;
            rho_r = pong[skip][IRHO];
            vx_r = pong[skip][IVX];
            vy_r = pong[skip][IVY];
            vz_r = pong[skip][IVZ];
            p_r = pong[skip][IP];
            a_r = pong[skip][IA];
            vn_r = vx_r * nx + vy_r * ny + vz_r * nz;
            if (moving_grid)
                gs = 0.5 * (nx * (wfl[pl][0] + wfl[pr][0]) +
                            ny * (wfl[pl][1] + wfl[pr][1]) +
                            nz * (wfl[pl][2] + wfl[pr][2]));
            vn_gs_l = vn_l - gs;
            vn_gs_r = vn_r - gs;
            var_l[IRHO] = 1;
            var_l[IRHOVX] = vx_l - vn_l * nx;
            var_l[IRHOVY] = vy_l - vn_l * ny;
            var_l[IRHOVZ] = vz_l - vn_l * nz;
            var_l[IRHOE] = ping[skip][IH];
            var_r[IRHO] = 1;
            var_r[IRHOVX] = vx_r - vn_r * nx;
            var_r[IRHOVY] = vy_r - vn_r * ny;
            var_r[IRHOVZ] = vz_r - vn_r * nz;
            var_r[IRHOE] = pong[skip][IH];
#define EXPAND_COUNT ZAP - 5
#include "expand.h"
            for (k = 5; k < ZAP; k++) {
                var_l[k] = ping[skip][k];
                var_r[k] = pong[skip][k];
            }
            ssw = (ping[skip][ISSW] + pong[skip][ISSW]) * 0.5;
            ssw_a = 1.0 - ssw;
            m_l = vn_gs_l / a_l;
            rhovn_l = p_12_l = 0.0;
            if (m_l >= 1) {
                rhovn_l = rho_l * vn_gs_l;
                p_12_l = p_l;
            } else if (m_l > -1) {
                const double mlp1 = 0.25 * SQR(m_l + 1);

                rhovn_l = rho_l * mlp1 * a_l;
                p_12_l = p_l * mlp1 * (2.0 - m_l);
            }
            m_r = vn_gs_r / a_r;
            rhovn_r = p_12_r = 0.0;
            if (m_r <= -1) {
                rhovn_r = rho_r * vn_gs_r;
                p_12_r = p_r;
            } else if (m_r < 1) {
                const double mrm1 = 0.25 * SQR(m_r - 1);

                rhovn_r = -rho_r * mrm1 * a_r;
                p_12_r = p_r * mrm1 * (2.0 + m_r);
            }
            p_12 = ssw * (p_12_l + p_12_r);
            rhovn_l *= ssw;
            rhovn_r *= ssw;
            rhovvn = rhovn_l * vn_l + rhovn_r * vn_r;
            a_m = MAX(a_l, a_r);
            if (fabs(vn_gs_l) < a_m) {
                double a_pl, m = vn_gs_l / a_m;

                p_pl = 0.25 * (2.0 - m) * SQR(m + 1) * p_l;
                a_pl =
                    0.5 * rho_r * p_l / ((rho_r * p_l + rho_l * p_r) *
                                         a_m);
                if (vn_gs_l >= 0)
                    vn_pl = vn_gs_l + a_pl * SQR(vn_gs_l - a_m);
                else
                    vn_pl = a_pl * SQR(vn_gs_l + a_m);
            } else {
                if (vn_gs_l >= 0.0) {
                    vn_pl = vn_gs_l;
                    p_pl = p_l;
                } else {
                    vn_pl = 0.0;
                    p_pl = 0.0;
                }
            }
            if (fabs(vn_gs_r) < a_m) {
                double a_mi, m = vn_gs_r / a_m;

                p_mi = 0.25 * (2.0 + m) * SQR(m - 1) * p_r;
                a_mi =
                    0.5 * rho_l * p_r / ((rho_r * p_l + rho_l * p_r) *
                                         a_m);
                if (vn_gs_r >= 0)
                    vn_mi = -a_mi * SQR(vn_gs_r - a_m);
                else
                    vn_mi = vn_gs_r - a_mi * SQR(vn_gs_r + a_m);
            } else {
                if (vn_gs_r >= 0.0) {
                    vn_mi = 0.0;
                    p_mi = 0.0;
                } else {
                    vn_mi = vn_gs_r;
                    p_mi = p_r;
                }
            }
            p_12 += (p_pl + p_mi) * ssw_a;
            vn_pl *= ssw_a;
            vn_mi *= ssw_a;
            rhovn = vn_pl * rho_l + vn_mi * rho_r;
            rhovn_la = rhovn_ra = 0.0;
            if (rhovn >= 0.0)
                rhovn_la = rhovn;
            else
                rhovn_ra = rhovn;
            rhovn_l += rhovn_la;
            rhovn_r += rhovn_ra;
            tmp = MIN(p_l, p_r);
            tmp = 10. * fabs(p_r - p_l) / tmp;
            if (1 <= tmp)
                s = 0.5;
            else
                s = 0.5 * tmp;
            rhovvn += (rhovn_la * vn_l + rhovn_ra * vn_r) * (0.5 - s);
            rhovvn +=
                (0.5 + s) * (rho_l * vn_l * vn_pl + rho_r * vn_r * vn_mi);
            atmp = ((vn_gs_l - a_l) < 0) && ((vn_gs_r - a_r) > 0);
            btmp = ((vn_gs_l + a_l) < 0) && ((vn_gs_r + a_r) > 0);
            if ((atmp && !btmp) || (!atmp && btmp)) {
                tmp = 0.125 * ssw_a;
                if (atmp && !btmp) {
                    tmp *= ((vn_gs_l - a_l) - (vn_gs_r - a_r));
                    rhovn_l -= tmp * rho_l;
                    rhovn_r += tmp * rho_r;
                } else {
                    tmp *= ((vn_gs_l + a_l) - (vn_gs_r + a_r));
                    rhovn_l -= tmp * rho_l;
                    rhovn_r += tmp * rho_r;
                }
            }
            rhovn_l *= area;
            rhovn_r *= area;
#define EXPAND_COUNT ZAP
#include "expand.h"
            for (k = 0; k < ZAP; k++)
                zap[k] = rhovn_l * var_l[k] + rhovn_r * var_r[k];
            rhovvn += p_12;
            zap[IRHOVX] += rhovvn * anx;
            zap[IRHOVY] += rhovvn * any;
            zap[IRHOVZ] += rhovvn * anz;
            zap[IRHOE] += gs * p_12 * area;
#define EXPAND_COUNT ZAP
#include "expand.h"
            for (k = 0; k < ZAP; k++) {
                dummy_4[pl][k] += zap[k];
                dummy_4[pr][k] -= zap[k];
            }
        }
    }
    free(ping);
    free(pong);
}

static void
kernel_2_2(DualGrid * grid, double dummy_1[][PR], double dummy_3[][AD],
           double dummy_2[][GR][3], double dummy_4[][CO])
{
    RangeList *cl;
    int (*hup)[2] = grid->hup;
    double (*hop)[3] = grid->hop,(*xx)[3] = grid->xx;
    const double two3 = 2. / 3.;
    const double junk_1 = 1.0;
    const double junk_2 = 1.0;
    const double junk_3 = 1.0;
    const double cgas = 1.0;

    for (cl = grid->fcl; cl != NULL; cl = cl->succ) {
        const int start = cl->start, stop = cl->stop, type = cl->type;
        int skip;

        if (type != TWO_SKIP) {
#include "nodep.h"
            for (skip = start; skip < stop; skip++) {
                double drho_dx, drho_dy, drho_dz, dvx_dx, dvx_dy, dvx_dz,
                    dvy_dx, dvy_dy, dvy_dz, dvz_dx, dvz_dy, dvz_dz, dp_dx,
                    dp_dy, dp_dz, dt_dx, dt_dy, dt_dz, tip_xx, tip_yy,
                    tip_zz, tip_xy, tip_xz, tip_yz, f, fold, fnew, dx, dy,
                    dz, rld;
                double zap[ZAP];
                const int p0 = hup[skip][0], p1 = hup[skip][1];
                const double nx = hop[skip][0], ny = hop[skip][1], nz =
                    hop[skip][2];
                const double rho0 = dummy_1[p0][IRHO], rho1 =
                    dummy_1[p1][IRHO], vx0 = dummy_1[p0][IVX], vx1 =
                    dummy_1[p1][IVX], vy0 = dummy_1[p0][IVY], vy1 =
                    dummy_1[p1][IVY], vz0 = dummy_1[p0][IVZ], vz1 =
                    dummy_1[p1][IVZ], pr0 = dummy_1[p0][IP], pr1 =
                    dummy_1[p1][IP], temp0 =
                    BLA_THREE(dummy_1, p0), temp1 = BLA_THREE(dummy_1, p1);
                const double rho = 0.5 * (rho0 + rho1), vx =
                    0.5 * (vx0 + vx1), vy = 0.5 * (vy0 + vy1), vz =
                    0.5 * (vz0 + vz1), p = 0.5 * (pr0 + pr1), rho_cg =
                    1 / (SQR(rho) * cgas);
                const double mue0 = BLA_ONE(temp0, dummy_1, p0), mue1 =
                    BLA_ONE(temp1, dummy_1, p1), kappa0 =
                    BLA_FOUR(temp0, dummy_1, p0), kappa1 =
                    BLA_FOUR(temp1, dummy_1, p1);
                const double mue_lam = 0.5 * (mue0 + mue1), kap_lam =
                    0.5 * (kappa0 + kappa1);
                double mue_eff = mue_lam, kappa = kap_lam, lambda;
                const double eddyv =
                    0.5 * (dummy_3[p0][AEDV] + dummy_3[p1][AEDV]);
                mue_eff += eddyv;
                kappa += grid->ptl_ratio * kap_lam / mue_lam * eddyv;
                drho_dx =
                    0.5 * (dummy_2[p0][IRHO][0] + dummy_2[p1][IRHO][0]);
                drho_dy =
                    0.5 * (dummy_2[p0][IRHO][1] + dummy_2[p1][IRHO][1]);
                drho_dz =
                    0.5 * (dummy_2[p0][IRHO][2] + dummy_2[p1][IRHO][2]);
                dvx_dx = 0.5 * (dummy_2[p0][IVX][0] + dummy_2[p1][IVX][0]);
                dvx_dy = 0.5 * (dummy_2[p0][IVX][1] + dummy_2[p1][IVX][1]);
                dvx_dz = 0.5 * (dummy_2[p0][IVX][2] + dummy_2[p1][IVX][2]);
                dvy_dx = 0.5 * (dummy_2[p0][IVY][0] + dummy_2[p1][IVY][0]);
                dvy_dy = 0.5 * (dummy_2[p0][IVY][1] + dummy_2[p1][IVY][1]);
                dvy_dz = 0.5 * (dummy_2[p0][IVY][2] + dummy_2[p1][IVY][2]);
                dvz_dx = 0.5 * (dummy_2[p0][IVZ][0] + dummy_2[p1][IVZ][0]);
                dvz_dy = 0.5 * (dummy_2[p0][IVZ][1] + dummy_2[p1][IVZ][1]);
                dvz_dz = 0.5 * (dummy_2[p0][IVZ][2] + dummy_2[p1][IVZ][2]);
                dp_dx = 0.5 * (dummy_2[p0][IP][0] + dummy_2[p1][IP][0]);
                dp_dy = 0.5 * (dummy_2[p0][IP][1] + dummy_2[p1][IP][1]);
                dp_dz = 0.5 * (dummy_2[p0][IP][2] + dummy_2[p1][IP][2]);
                dt_dx = (dp_dx * rho - drho_dx * p) * rho_cg;
                dt_dy = (dp_dy * rho - drho_dy * p) * rho_cg;
                dt_dz = (dp_dz * rho - drho_dz * p) * rho_cg;
                if (grid->level == LOW_LEVEL) {
                    dx = xx[p1][0] - xx[p0][0];
                    dy = xx[p1][1] - xx[p0][1];
                    dz = xx[p1][2] - xx[p0][2];
                    rld = 1 / sqrt(dx * dx + dy * dy + dz * dz);
                    dx *= rld;
                    dy *= rld;
                    dz *= rld;
                    fold = dvx_dx * dx + dvx_dy * dy + dvx_dz * dz;
                    fnew = (vx1 - vx0) * rld;
                    f = fnew - fold;
                    dvx_dx += f * dx;
                    dvx_dy += f * dy;
                    dvx_dz += f * dz;
                    fold = dvy_dx * dx + dvy_dy * dy + dvy_dz * dz;
                    fnew = (vy1 - vy0) * rld;
                    f = fnew - fold;
                    dvy_dx += f * dx;
                    dvy_dy += f * dy;
                    dvy_dz += f * dz;
                    fold = dvz_dx * dx + dvz_dy * dy + dvz_dz * dz;
                    fnew = (vz1 - vz0) * rld;
                    f = fnew - fold;
                    dvz_dx += f * dx;
                    dvz_dy += f * dy;
                    dvz_dz += f * dz;
                    fold = dt_dx * dx + dt_dy * dy + dt_dz * dz;
                    fnew = (temp1 - temp0) * rld;
                    f = fnew - fold;
                    dt_dx += f * dx;
                    dt_dy += f * dy;
                    dt_dz += f * dz;
                }
                lambda = -two3 * mue_eff;
                tip_xx = lambda * (dvy_dy + dvz_dz - 2.0 * dvx_dx);
                tip_yy = lambda * (dvx_dx + dvz_dz - 2.0 * dvy_dy);
                tip_zz = lambda * (dvx_dx + dvy_dy - 2.0 * dvz_dz);
                tip_xy = mue_eff * (dvx_dy + dvy_dx);
                tip_xz = mue_eff * (dvx_dz + dvz_dx);
                tip_yz = mue_eff * (dvy_dz + dvz_dy);
                zap[IVX] = -(tip_xx * nx + tip_xy * ny + tip_xz * nz);
                zap[IVY] = -(tip_xy * nx + tip_yy * ny + tip_yz * nz);
                zap[IVZ] = -(tip_xz * nx + tip_yz * ny + tip_zz * nz);
                zap[IRHOE] =
                    -((vx * tip_xx + vy * tip_xy + vz * tip_xz +
                       kappa * dt_dx) * nx + (vx * tip_xy + vy * tip_yy +
                                              vz * tip_yz +
                                              kappa * dt_dy) * ny +
                      (vx * tip_xz + vy * tip_yz + vz * tip_zz +
                       kappa * dt_dz) * nz);
                dummy_4[p0][IRHOVX] += zap[IVX];
                dummy_4[p0][IRHOVY] += zap[IVY];
                dummy_4[p0][IRHOVZ] += zap[IVZ];
                dummy_4[p0][IRHOE] += zap[IRHOE];
                dummy_4[p1][IRHOVX] -= zap[IVX];
                dummy_4[p1][IRHOVY] -= zap[IVY];
                dummy_4[p1][IRHOVZ] -= zap[IVZ];
                dummy_4[p1][IRHOE] -= zap[IRHOE];
            }
        }
    }
}

static void
kernel_2_3(DualGrid * grid, double dummy_1[][PR], double dummy_3[][AD],
           double dummy_2[][GR][3], double dummy_4[][CO])
{
    int pnt;
    double source;
    const int npnt = grid->nallpoints;
    const double *pvolume = grid->pvolume;
    const double *wdist = grid->wdist;
    const double rscale = 1. / tanh(1.0);
    const double four3 = 4. / 3.;
    const double cv1p3 = 1.0;
    const double rkap2 = 1.0;
    const double cw2 = 1.0;
    const double cw3p6 = 1.0;
    const double cb1 = 1.0;
    const double cb2s = 1.0;
    const double cw1k = 1.0;
    const double junk_1 = 1.0;
    const double junk_3 = 1.0;

#include "nodep.h"
    for (pnt = 0; pnt < npnt; pnt++) {
        double diffsource = 0;
        double production = 0;
        double wall_dest = 0;
        const double rho = dummy_1[pnt][IRHO];
        const double nuet = dummy_1[pnt][INUE];
        const double muet = rho * nuet;
        const double temp = BLA_THREE(dummy_1, pnt);
        const double muel = BLA_ONE(temp, dummy_1, pnt);
        const double xsi = muet / muel;
        const double xsi3 = xsi * xsi * xsi;
        const double fv1 = xsi3 / (xsi3 + cv1p3);
        const double s_fac = fv1 + 1 / MAX(xsi, EPS);
        const double dist = wdist[pnt];
        const double d2 = dist * dist;
        const double dudx = dummy_2[pnt][IVX][0];
        const double dudy = dummy_2[pnt][IVX][1];
        const double dudz = dummy_2[pnt][IVX][2];
        const double dvdx = dummy_2[pnt][IVY][0];
        const double dvdy = dummy_2[pnt][IVY][1];
        const double dvdz = dummy_2[pnt][IVY][2];
        const double dwdx = dummy_2[pnt][IVZ][0];
        const double dwdy = dummy_2[pnt][IVZ][1];
        const double dwdz = dummy_2[pnt][IVZ][2];
        const double tmp =
            (SQR(dudx) + SQR(dvdy) + SQR(dwdz)) * four3 + SQR(dudy) +
            SQR(dudz) + SQR(dvdx) + SQR(dvdz) + SQR(dwdx) + SQR(dwdy) +
            dudy * dvdx + dudz * dwdx + dvdz * dwdy;
        const double strain = sqrt(MAX(tmp, 0.));
        const double dnue_dx = dummy_2[pnt][INUE][0];
        const double dnue_dy = dummy_2[pnt][INUE][1];
        const double dnue_dz = dummy_2[pnt][INUE][2];
        const double gradnue2 =
            dnue_dx * dnue_dx + dnue_dy * dnue_dy + dnue_dz * dnue_dz;
        const double rs = nuet * rkap2 / d2;
        const double muet_rs = rho * nuet * rs;
        double fw, r, g, s, lim_g;

        s = MAX(strain * s_fac, EPS);
        r = tanh(rs / s) * rscale;
        g = r + cw2 * (P6(r) - r);
        lim_g = (1 + cw3p6) / (P6(g) + cw3p6);
        fw = g * SQRT_SIX(lim_g);
        production = cb1 * strain * s_fac * muet;
        diffsource = cb2s * rho * gradnue2;
        wall_dest = cw1k * muet_rs * fw;
        source = production - wall_dest;
        source += diffsource;
        dummy_4[pnt][IRHONUE] -= source * pvolume[pnt];
}} static void

kernel_2_4(DualGrid * grid, double dummy_1[][PR], double dummy_3[][AD],
           double dummy_2[][GR][3], double dummy_4[][CO])
{
    int (*hup)[2] = grid->hup;
    double (*hop)[3] = grid->hop,(*xx)[3] = grid->xx;
    RangeList *cl;
    const double rsigm = 1.0;
    const double junk_1 = 1.0;
    const double junk_3 = 1.0;

    for (cl = grid->fcl; cl != NULL; cl = cl->succ) {
        const int start = cl->start, stop = cl->stop, type = cl->type;
        int skip;

        if (type != TWO_SKIP) {
#include "nodep.h"
            for (skip = start; skip < stop; skip++) {
                double dnue_dx, dnue_dy, dnue_dz, f, fold, fnew, dx, dy,
                    dz, rld, zap;
                const int p0 = hup[skip][0], p1 = hup[skip][1];
                const double nx = hop[skip][0], ny = hop[skip][1], nz =
                    hop[skip][2], rho0 = dummy_1[p0][IRHO], rho1 =
                    dummy_1[p1][IRHO], ed0 = dummy_1[p0][INUE], ed1 =
                    dummy_1[p1][INUE];
                const double mue_t = 0.5 * (rho0 * ed0 + rho1 * ed1);
                const double tem0 = BLA_THREE(dummy_1, p0), tem1 =
                    BLA_THREE(dummy_1, p1);
                const double mue_l0 = BLA_ONE(tem0, dummy_1, p0), mue_l1 =
                    BLA_ONE(tem1, dummy_1, p1), mue_lam =
                    0.5 * (mue_l0 + mue_l1), mue_eff =
                    (mue_lam + mue_t) * rsigm;
                dnue_dx =
                    0.5 * (dummy_2[p0][INUE][0] + dummy_2[p1][INUE][0]);
                dnue_dy =
                    0.5 * (dummy_2[p0][INUE][1] + dummy_2[p1][INUE][1]);
                dnue_dz =
                    0.5 * (dummy_2[p0][INUE][2] + dummy_2[p1][INUE][2]);
                if (grid->level == LOW_LEVEL) {
                    dx = xx[p1][0] - xx[p0][0];
                    dy = xx[p1][1] - xx[p0][1];
                    dz = xx[p1][2] - xx[p0][2];
                    rld = 1 / sqrt(dx * dx + dy * dy + dz * dz);
                    dx *= rld;
                    dy *= rld;
                    dz *= rld;
                    fold = dnue_dx * dx + dnue_dy * dy + dnue_dz * dz;
                    fnew = (ed1 - ed0) * rld;
                    f = fnew - fold;
                    dnue_dx += f * dx;
                    dnue_dy += f * dy;
                    dnue_dz += f * dz;
                }
                zap =
                    -(dnue_dx * nx + dnue_dy * ny +
                      dnue_dz * nz) * mue_eff;
                dummy_4[p0][IRHONUE] += zap;
                dummy_4[p1][IRHONUE] -= zap;
            }
        }
    }
}
