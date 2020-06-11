#ifndef UTIL_H
#define UTIL_H

#include <math.h>

#ifndef SQR
#define SQR(x) ((x) * (x))
#endif

#ifndef MIN
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#endif

#ifndef MAX
#define MAX(a, b) ( (a) > (b) ? (a) : (b))
#endif

#ifndef MAX0
#define MAX0(a) ((a) > 0 ? (a) : 0)
#endif

#ifndef MIN0
#define MIN0(a) ((a) < 0 ? (a) : 0)
#endif

#ifndef ABS
#define ABS(a) ((a) < 0 ? -1 * (a) : (a) )
#endif

#ifndef ISIGN
#define ISIGN(a) ((a) < 0 ? -1  : 1 )
#endif

#ifndef POW
#define POW(x, y) exp( (y) * log(x) )
#endif


#ifndef POW_LIM
#define POW_LIM(x, y) exp( (y) * log(MAX((x), 1e-30)) )
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


#define EPS             1.e-16
#define P6(x)           ((x) * (x) * (x) * (x) * (x) * (x))
#define SQRT_SIX(x)     (exp(0.166666666666666 * log((x))))
#define LIM_G_MIN       1.e-32
#define OMEGA 1.0
#define POW_ZETA(x) sqrt(x)
#define FOUR_THIRDS (4.0 / 3.0)


#define BLA_ONE(temp, dummy_1, pnt)\
         (junk_1 * BLA_FIVE(temp))

#define BLA_TWO(hop, skip)  \
        (sqrt(SQR(hop[skip][0])+SQR(hop[skip][1])+SQR(hop[skip][2])))

#define BLA_THREE(dummy_1, pnt) \
         (dummy_1[pnt][IP] / dummy_1[pnt][IRHO])

#define BLA_FOUR(temp, dummy_1, pnt)\
         (junk_2 * BLA_FIVE(temp))

#define BLA_FIVE(temp) \
         ((temp) * sqrt((temp)) * (1.0 + junk_3)\
                                / ((temp) + junk_3))
#define BLA_SIX(dummy_1, pnt) \
         (junk_4)

#define BLA_SEVEN(dummy_1, pnt) \
         (sqrt(junk_5 * dummy_1[pnt][IP] / dummy_1[pnt][IRHO]))


#endif
