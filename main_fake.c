
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <string.h>
#include <sys/resource.h>
#include <sys/time.h>
#include "types.h"
#include "lim_fake.h"
#include "flux_fake.h"
#include "smooth_fake.h"
#include "util.h"
#include "flop.h"
double second();
void *test_malloc(int pnt);
static void usage(int id);
static void shuffle(int size, int np, int id, double *in, double *out);
static void free_struc(DualGrid * grid, WorkSpace * work);
static void alloc(DualGrid * grid, WorkSpace * work, double dummy_1[][PR],
                  double dummy_2[][GR][3], double dummy_3[][AD],
                  double dummy_4[][CO], int npoints, int nskips, int ncs,
                  int id);
static void init(DualGrid * grid, WorkSpace * work, double dummy_1[][PR],
                 double dummy_2[][GR][3], double dummy_3[][AD],
                 double dummy_4[][CO], int npoints, int nskips, int ncs);
#define FW 4.0
#define CW 4000.0
#define RC 12
#define LI_ST 2
#define ZA_ST 1
#define SM_ST 2
#define SC 1.0
#define SCP 1.0
#define SM 2
#define RE 2
#define IDE 16
int
main(int argc, char **argv)
{
    DualGrid *grid = NULL;
    WorkSpace *work = NULL;
    double (*dummy_1)[PR] = NULL;
    double (*dummy_2)[GR][3] = NULL;
    double (*dummy_3)[AD] = NULL;
    double (*dummy_4)[CO] = NULL;
    double *out = NULL, *in = NULL;
    double a1 = 0.0, a2 = 0.0, a3 = 0.0;
    int np = 0, id = 0, pnts = 262144, i, j;
    int npoints = 0, nskips = 0, ncs = 0;
    int size = 0;
    int steps = 10;
    int prn = 0;
    static double sa1 = 0.0, sa2 = 0.0;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    if (argc > 1) {
        if (argc != 5)
            usage(id);
        if (strcmp(argv[1], "-n") || strcmp(argv[3], "-s"))
            usage(id);
        if (!strcmp(argv[1], "-n")) {
            pnts = atoi(argv[2]);
        }
        if (!strcmp(argv[3], "-s")) {
            steps = atoi(argv[4]);
        }
    }
    npoints = pnts;
    nskips = (int) npoints *FW;

    ncs = (int) nskips / CW;
    grid = (DualGrid *) test_malloc(sizeof(DualGrid));;
    work = (WorkSpace *) test_malloc(sizeof(WorkSpace));;
    dummy_1 = (double (*)[PR]) test_malloc(npoints * sizeof(double[PR]));
    dummy_2 =
        (double (*)[GR][3]) test_malloc(npoints * sizeof(double[GR][3]));
    dummy_3 = (double (*)[AD]) test_malloc(npoints * sizeof(double[AD]));
    dummy_4 = (double (*)[CO]) test_malloc(npoints * sizeof(double[CO]));
    alloc(grid, work, dummy_1, dummy_2, dummy_3, dummy_4, npoints, nskips,
          ncs, id);
    size = (int) pow((double) npoints, 0.66666666) * 444.0 / RC;
    out = (double *) test_malloc(size * sizeof(double));
    in = (double *) test_malloc(size * sizeof(double));
    if (id == 0) {
        printf("This is TauBench.\n");
        printf("Evaluating kernels - please be patient.\n");
    }
    for (i = 0; i < steps; i++) {
        if (id == 0) {
            printf(".");
            if (i == steps - 1) {
                printf("\n\n");
            }
            fflush(stdout);
            a1 = second();
        }
        init(grid, work, dummy_1, dummy_2, dummy_3, dummy_4, npoints,
             nskips, ncs);
        for (j = 0; j < LI_ST; j++) {
            const double sc = SC;
            const double scp = SCP;

            if (i == steps - 1 && j == LI_ST - 1) {
                prn = 1;
            } else {
                prn = 0;
            }
            kernel_1_0(grid, work, dummy_1, dummy_2, sc, scp, id, prn);
        }
        for (j = 0; j < ZA_ST; j++) {
            if (i == steps - 1 && j == ZA_ST - 1) {
                prn = 1;
            } else {
                prn = 0;
            }
            kernel_2_0(grid, work, dummy_1, dummy_3, dummy_2, dummy_4, id,
                       prn);
        }
        for (j = 0; j < SM_ST; j++) {
            int sm = SM;
            int re = RE;

            if (i == steps - 1 && j == SM_ST - 1) {
                prn = 1;
            } else {
                prn = 0;
            }
            kernel_3_0(grid, work, dummy_4, dummy_3, dummy_1, sm, re, id,
                       prn);
        }
        if (id == 0) {
            a2 = second();
            sa1 += (a2 - a1);
        }
        shuffle(size, np, id, in, out);
        if (id == 0) {
            a3 = second();
            sa2 += (a3 - a2);
        }
    }
    if (id == 0) {
        printf("\n");
        printf("               total : %10.3f secs - %10.3f mflops\n\n",
               sa1 + sa2, np * TZE * npoints / (sa1 + sa2) * steps);
        printf("points     : %10d\n", pnts);
        printf("steps      : %10d\n", steps);
        printf("procs      : %10d\n\n", np);
        printf("comp       : %10.3f secs\n", sa1);
        printf("comm       : %10.3f secs\n", sa2);
        printf("comm ratio : %10.3f\n\n", sa2 / sa1);
    }
    free_struc(grid, work);
    free(grid);
    free(work);
    free(in);
    free(out);
    free(dummy_1);
    free(dummy_2);
    free(dummy_3);
    free(dummy_4);
    MPI_Finalize();
    return 0;
}

static void
alloc(DualGrid * grid, WorkSpace * work, double dummy_1[][PR],
      double dummy_2[][GR][3], double dummy_3[][AD], double dummy_4[][CO],
      int npoints, int nskips, int ncs, int id)
{
    int (*hup)[2] = NULL;
    int cl = 0;
    int kl = 0;
    RangeList *init = NULL, **succ = NULL;

    work->wfl = (double (*)[3]) test_malloc(npoints * sizeof(double[3]));
    work->plim =
        (double (*)[LIM]) test_malloc(npoints * sizeof(double[LIM]));
    grid->xx = (double (*)[3]) test_malloc(npoints * sizeof(double[3]));
    grid->pvolume = (double *) test_malloc(npoints * sizeof(double));
    grid->wdist = (double *) test_malloc(npoints * sizeof(double));
    grid->hup = (int (*)[2]) test_malloc(nskips * sizeof(int[2]));
    grid->hop = (double (*)[3]) test_malloc(nskips * sizeof(double[3]));
    grid->hip = (double *) test_malloc(nskips * sizeof(double));
    grid->fcl = NULL;
    succ = &grid->fcl;
    for (cl = 0; cl <= ncs; cl++) {
        init = (RangeList *) test_malloc(sizeof(RangeList));
        init->start = kl;
        kl += (int) CW;
        kl = MIN(nskips, kl);
        init->stop = kl;
        init->type = 0;
        init->succ = NULL;
        *succ = init;
        succ = &init->succ;
    } hup = grid->hup;
    for (init = grid->fcl; init != NULL; init = init->succ) {
        const int start = init->start, stop = init->stop;
        int skip, k = 0, ili = 0, iri = 0, iki = 1, idi = 1;
        int isi, ifi, c1, c2;

        isi = (int) start / FW;
        ifi = (int) stop / FW;
        for (skip = start; skip < stop; skip++) {
            if (k == (int) FW) {
                k = 0;
                ili++;
                idi++;
                if (idi > IDE) {
                    idi = 1;
                }
            }
            k++;
            iri = ili + iki;
            c1 = isi + ili;
            if (c1 >= ifi) {
                c1 = isi + (ili % (ifi - isi));
            }
            c2 = isi + iri;
            if (c2 >= ifi) {
                c2 = isi + (iri % (ifi - isi));
            }
            hup[skip][0] = c1;
            hup[skip][1] = c2;
            iki += idi;
        }
    }
}

static void
init(DualGrid * grid, WorkSpace * work, double dummy_1[][PR],
     double dummy_2[][GR][3], double dummy_3[][AD], double dummy_4[][CO],
     int npoints, int nskips, int ncs)
{
    double (*hop)[3] = NULL;
    double *hip = NULL;
    double (*wfl)[3] = NULL;
    double (*plim)[LIM] = NULL;
    double (*xx)[3] = NULL;
    double *pvolume = NULL;
    double *wdist = NULL;
    int pnt, skip;
    int i, j;

    grid->nallpoints = npoints;
    grid->nallskips = nskips;
    grid->level = LOW_LEVEL;
    grid->moving = 0;
    grid->ptl_k = 1.0;
    grid->ptl_ratio = 1.0;
    wfl = work->wfl;
    plim = work->plim;
    xx = grid->xx;
    pvolume = grid->pvolume;
    wdist = grid->wdist;
    hop = grid->hop;
    hip = grid->hip;
    for (pnt = 0; pnt < npoints; pnt++) {
        for (i = 0; i < 3; i++) {
            wfl[pnt][i] = (double) (i + 1);
    }} for (pnt = 0; pnt < npoints; pnt++) {
        for (i = 0; i < LIM; i++) {
            plim[pnt][i] = (double) (i + 1);
    }} for (pnt = 0; pnt < npoints; pnt++) {
        for (i = 0; i < 3; i++) {
            xx[pnt][i] = (double) (i + 1);
    }} for (pnt = 0; pnt < npoints; pnt++) {
        pvolume[pnt] = 1.0;
    }
    for (pnt = 0; pnt < npoints; pnt++) {
        wdist[pnt] = 1.0;
    }
    for (pnt = 0; pnt < npoints; pnt++) {
        for (i = 0; i < PR; i++) {
            dummy_1[pnt][i] = (double) (i + 1);
    }} for (pnt = 0; pnt < npoints; pnt++) {
        for (i = 0; i < GR; i++) {
            for (j = 0; j < 3; j++) {
                dummy_2[pnt][i][j] = (double) (j + 1);
    }}} for (pnt = 0; pnt < npoints; pnt++) {
        for (i = 0; i < AD; i++) {
            dummy_3[pnt][i] = 1.0;
        }
    }
    for (pnt = 0; pnt < npoints; pnt++) {
        for (i = 0; i < CO; i++) {
            dummy_4[pnt][i] = 1.0;
        }
    }
    for (skip = 0; skip < nskips; skip++) {
        for (i = 0; i < 3; i++) {
            hop[skip][i] = (double) (i + 1);
    }} for (skip = 0; skip < nskips; skip++) {
        hip[skip] = 1.0;
    }
}

static void
usage(int id)
{
    if (id == 0) {
        fprintf(stderr,
                "usage: taubench -n <gridpoints> -s <pseudosteps>\n");
    }
    exit(1);
}

static void
shuffle(int size, int np, int id, double *in, double *out)
{
    int dest, src, i;
    double ssum = 1.0;
    double rsum = 1.0;
    MPI_Status status[2 * RC];
    MPI_Request msgid[2 * RC];

    for (i = 0; i < size; i++) {
        out[i] = 1.0;
    }
    for (i = 0; i < RC; i++) {
        MPI_Allreduce(&ssum, &rsum, 1, MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD);
        dest = (id == np - 1) ? 0 : id + 1;
        MPI_Isend(out, size, MPI_DOUBLE, dest, 5, MPI_COMM_WORLD,
                  &msgid[i * 2 + 0]);
        src = (id == 0) ? np - 1 : id - 1;
        MPI_Irecv(in, size, MPI_DOUBLE, src, 5, MPI_COMM_WORLD,
                  &msgid[i * 2 + 1]);
    }
    MPI_Waitall(2 * RC, msgid, status);
}

static void
free_struc(DualGrid * grid, WorkSpace * work)
{
    int id;
    RangeList *init = NULL, *succ = NULL;

    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    free(work->wfl);
    free(work->plim);
    free(grid->xx);
    free(grid->pvolume);
    free(grid->wdist);
    free(grid->hop);
    free(grid->hip);
    free(grid->hup);
    init = grid->fcl;
    while (init != NULL) {
        succ = init->succ;
        free(init);
        init = succ;
    }
}

void *
test_malloc(int pnt)
{
    void *try = NULL;

    try = malloc(pnt);
    if (try == NULL) {
        printf("out of memory.\n");
        exit(1);
    }
    return try;
}

double
second()
{
    double tarray[3];
    double t1;
    struct rusage RU;
    struct timeval tv;

    getrusage(RUSAGE_SELF, &RU);
    gettimeofday(&tv, (struct timezone *) 0);
    tarray[0] = RU.ru_utime.tv_sec + (double) RU.ru_utime.tv_usec * 1e-6;
    t1 = tarray[0];
    return t1;
}
