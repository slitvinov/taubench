#include "types.h"
#include "lim.h"
#include "flop.h"
#include "main.h"
#include "util.h"
#include <stdio.h>
#include <stdlib.h>
static void kernel_1_1(DualGrid *grid, double dummy_1[][PR],
                       double dumin[][LIM], double dumax[][LIM], double eps,
                       int id);
void kernel_1_0(DualGrid *grid, WorkSpace *work, double dummy_1[][PR],
                double dummy_2[][GR][3], double sc, double scp, int id,
                int prn) {
  const double eps = 1.0e-08;
  const double threshold = 1.0 - eps;
  const int npnt = grid->nallpoints;
  double(*xx)[3] = grid->xx;
  double(*plim)[LIM] = work->plim;
  double(*dumax)[LIM], (*dumin)[LIM], (*newlim)[LIM];
  int(*hup)[2] = grid->hup;
  int pnt, skip;
  int k;
  static double st1 = 0.0, st2 = 0.0;
  static int stc = 0;
  double t1 = 0.0, t2 = 0.0, t3 = 0.0, t4 = 0.0;
  RangeList *cl;

  if (id == 0) {
    t1 = second();
  }
  dumin = (double(*)[LIM])test_malloc(npnt * sizeof(double[LIM]));
  ;
  dumax = (double(*)[LIM])test_malloc(npnt * sizeof(double[LIM]));
  ;
  newlim = (double(*)[LIM])test_malloc(npnt * sizeof(double[LIM]));
  ;
  if (id == 0) {
    t3 = second();
  }
  kernel_1_1(grid, dummy_1, dumin, dumax, eps, id);
  if (id == 0) {
    t4 = second();
  }
#include "nodep.h"
  for (pnt = 0; pnt < npnt; pnt++) {
#define EXPAND_COUNT LIM
#include "expand.h"
    for (k = 0; k < LIM; k++) {
      newlim[pnt][k] = 2;
    }
  }
  for (cl = grid->fcl; cl != NULL; cl = cl->succ) {
    const int start = cl->start;
    const int stop = cl->stop;
    const int type = cl->type;

    if (type != TWO_SKIP) {
#include "nodep.h"
      for (skip = start; skip < stop; skip++) {
        const int p0 = hup[skip][0], p1 = hup[skip][1];
        const double dx = 0.5 * (xx[p1][0] - xx[p0][0]);
        const double dy = 0.5 * (xx[p1][1] - xx[p0][1]);
        const double dz = 0.5 * (xx[p1][2] - xx[p0][2]);

#define EXPAND_COUNT LIM
#include "expand.h"
        for (k = 0; k < LIM; k++) {
          double r0, r1;
          const double newlim0 = newlim[p0][k], newlim1 = newlim[p1][k];
          const double d0 = dx * dummy_2[p0][k][0] + dy * dummy_2[p0][k][1] +
                            dz * dummy_2[p0][k][2];
          const double d1 = -dx * dummy_2[p1][k][0] - dy * dummy_2[p1][k][1] -
                            dz * dummy_2[p1][k][2];
          if (d0 == 0) {
            r0 = 2;
          } else {
            const double rd0 = 1 / d0;

            r0 = (d0 > 0) ? (dumax[p0][k] * rd0) : (dumin[p0][k] * rd0);
          }
          if (d1 == 0) {
            r1 = 2;
          } else {
            const double rd1 = 1 / d1;

            r1 = (d1 > 0) ? (dumax[p1][k] * rd1) : (dumin[p1][k] * rd1);
          }
          newlim[p0][k] = MIN(newlim0, r0);
          newlim[p1][k] = MIN(newlim1, r1);
        }
      }
    }
  }
#include "nodep.h"
  for (pnt = 0; pnt < npnt; pnt++) {
#define EXPAND_COUNT LIM
#include "expand.h"
    for (k = 0; k < LIM; k++) {
      const double r = newlim[pnt][k], rr = r * r;
      const double tmplim = (rr + 2.0 * r) / (rr + r + 2.0);

      if (sc > threshold)
        plim[pnt][k] = tmplim;
      else {
        double d_lim = plim[pnt][k] - tmplim;

        if (d_lim < 0)
          plim[pnt][k] -= d_lim * scp;
        else
          plim[pnt][k] -= d_lim * sc;
      }
    }
  }
  free(newlim);
  free(dumax);
  free(dumin);
  if (id == 0) {
    t2 = second();
    st1 += (t2 - t1);
    st2 += (t4 - t3);
    stc++;
    if (prn != 0) {
      printf("        - kernel_1_0 : %10.3f secs - %10.3f mflops\n", st1,
             GZE * npnt / st1 * stc);
      printf("        - kernel_1_1 : %10.3f secs - %10.3f mflops\n", st2,
             GON * npnt / st2 * stc);
    }
  }
}

static void kernel_1_1(DualGrid *grid, double dummy_1[][PR],
                       double dumin[][LIM], double dumax[][LIM], double eps,
                       int id) {
  const int npnt = grid->nallpoints;
  int(*hup)[2] = grid->hup;
  int pnt, skip, k;
  RangeList *cl;

#include "nodep.h"
  for (pnt = 0; pnt < npnt; pnt++) {
#define EXPAND_COUNT LIM
#include "expand.h"
    for (k = 0; k < LIM; k++) {
      dumax[pnt][k] = eps;
      dumin[pnt][k] = -eps;
    }
  }
  for (cl = grid->fcl; cl != NULL; cl = cl->succ) {
    const int start = cl->start;
    const int stop = cl->stop;
    const int type = cl->type;

    if (type != TWO_SKIP) {
#include "nodep.h"
      for (skip = start; skip < stop; skip++) {
        const int p0 = hup[skip][0], p1 = hup[skip][1];

#define EXPAND_COUNT LIM
#include "expand.h"
        for (k = 0; k < LIM; k++) {
          const double du = dummy_1[p1][k] - dummy_1[p0][k];
          const double dum = -du;
          const double dumin0 = dumin[p0][k], dumax0 = dumax[p0][k];
          const double dumin1 = dumin[p1][k], dumax1 = dumax[p1][k];
          dumin[p0][k] = MIN(dumin0, du);
          dumax[p0][k] = MAX(dumax0, du);
          dumin[p1][k] = MIN(dumin1, dum);
          dumax[p1][k] = MAX(dumax1, dum);
        }
      }
    }
  }
}
