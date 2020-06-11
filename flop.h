#ifdef ITASCALE
/* Intel Itanium Scale */
#define GZE 0.0006467742
#define GON 0.0000886416
#define ZON 0.0003328242
#define ZTW 0.0015767801
#define ZTH 0.0003937917
#define ZFO 0.0008047718
#define SAL 0.0027712333
#define TZE 0.0080839163
#define _SCALE_
#endif

#ifdef SR8KSCALE
/* Hitachi SR8000 Scale */
#define GZE 0.0005044010
#define GON 0.0002004500
#define ZON 0.0006857705
#define ZTW 0.0008298101
#define ZTH 0.0002496210
#define ZFO 0.0002405571
#define SAL 0.0004410886
#define TZE 0.0039511980
#define _SCALE_
#endif

#ifndef _SCALE_
#define SX6SCALE
#endif

#ifdef SX6SCALE
/* Nec SX6 Scale */
#define GZE 0.0007035689
#define GON 0.0001090219
#define ZON 0.0009931064
#define ZTW 0.0012199851
#define ZTH 0.0002462927
#define ZFO 0.0005544819
#define SAL 0.0013693298
#define TZE 0.0063305670
#define _SCALE_
#endif
