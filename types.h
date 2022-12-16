#define CO 6
#define ZAP 6
#define PR 6
#define GR 6
#define ST 15
#define AD 1
#define LIM 5

typedef struct RangleList_Struct {

  int start, stop, type;
  int vector_coloring;

  struct RangleList_Struct *succ;
} RangeList;

typedef struct {
  double (*wfl)[3];
  double (*plim)[LIM];

} WorkSpace;

typedef struct {

  int nallpoints;
  double (*xx)[3];
  double *pvolume;
  double *wdist;
  int nallskips;
  int moving;
  int level;
  double ptl_k;
  double ptl_ratio;
  int (*hup)[2];
  int (*fngb)[2];
  double (*hop)[3];
  double *hip;

  RangeList *fcl;

} DualGrid;

#define TWO_SKIP 1
#define LOW_LEVEL 1

#define IRHO 0
#define IRHOVX 1
#define IRHOVY 2
#define IRHOVZ 3
#define IRHOE 4
#define IRHONUE 5
#define ISSW 6
#define IT 7
#define IA 8
#define IH 9
#define ICV 10
#define IMUE 11
#define IKAPPA 12
#define IS 13
#define ICP 14

#define AEDV 0

#define IVX IRHOVX
#define IVY IRHOVY
#define IVZ IRHOVZ
#define IP IRHOE
#define INUE IRHONUE
#define IPEX ISSW
#define ILIM ISSW
