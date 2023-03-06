#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef DNACONSTANTS_H_
#define DNACONSTANTS_H_

#define DNA_VERSION_NUM (0.31)
static const char DNA_RELEASE_DATE[] = "17 Feb 2023";

/**--------------------------------------------------------- 
Not all of the following hash functions for pre-defined constants and
operations are currently in use, but might be useful for later work.
The header <DNA-constants.h> must be included to access the hash functions.
---------------------------------------------------------**/

typedef double DNA_FLOAT;
#define DNA_EPS (1.0e-15)
#define DNA_SMALL (1.0e-30)
#define DNA_POW(a, b) (pow(a, b))
#define DNA_SQRT(a) (sqrt(a))
#define DNA_ABS(a) (fabs(a))
#define DNA_SIN(a) (sin(a))
#define DNA_COS(a) (cos(a))
#define DNA_EXP(a) (exp(a))
#define DNA_LOG(a) (log(a))
#define DNA_CEIL(a) (ceil(a))
#define DNA_STRINGTOFLOAT(a) (strtod(a, &dummy))

/** Simple operations **/
#define DNA_POW2(x) ((x) * (x))
#define DNA_POW3(x) ((x) * (x) * (x))
#define DNA_POW4(x) ((x) * (x) * (x) * (x))
#define DNA_MAX(a, b) ((a) > (b) ? (a) : (b))
#define DNA_MAX3(a, b, c) (((a) >= (b) && (a) >= (c)) ? (a) : ((b) >= (a) && (b) >= (c)) ? (b) : (c))
#define DNA_MIN(a, b) ((a) < (b) ? (a) : (b))
#define DNA_MIN3(a, b, c) (((a) <= (b) && (a) <= (c)) ? (a) : ((b) <= (a) && (b) <= (c)) ? (b) : (c))

/** Pre-defined constants **/
#define DNA_ONETHIRD (0.3333333333333333333333333333333333333333)
#define DNA_TWOTHIRD (0.6666666666666666666666666666666666666667)
#define DNA_PI (3.1415926535897932384626433832795028841972)
#define DNA_E (2.71828182845904523536028747135266249775724709369995)

/** I/O Priority (bit-wise comparison) **/
#define DNA_LOWPRIO (1)
#define DNA_MEDIUMPRIO (2)
#define DNA_HIGHPRIO (4)

/** Misc **/
#define DNA_DATA_ALLOC_INCREMENT (10000)
#define DNA_STRINGLENGTH (512)
#define DNA_STRINGLENGTH_LONG (1024)
#define DNA_STRINGLENGTH_SPRINTF (1024)
#define DNA_STRINGLENGTH_SPRINTF_LONG (2048)

/** Debugging aid **/
#define DNA_WHERE printf("HERE - %s:%d\n", __FILE__, __LINE__);
#define DNA_WHERE_INT(a)                                \
  printf("HERE - %s:%d = %i\n", __FILE__, __LINE__, a); \
  <
#endif
