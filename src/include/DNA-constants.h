#include <math.h>

#ifndef DNACONSTANTS_H_
#define DNACONSTANTS_H_

#define DNA_VERSION_NUM (1.0)
static const char DNA_RELEASE_DATE[] = "7 Apr 2023";

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
#define DNA_POW2(x) ((x) * (x))
#define DNA_POW3(x) ((x) * (x) * (x))
#define APECSS_MIN(a, b) ((a) < (b) ? (a) : (b))
#define DNA_PI (3.1415926535897932384626433832795028841972)
#define DNA_STRINGLENGTH (512)
#define DNA_STRINGLENGTH_SPRINTF (1024)

/** Debugging aid **/
#define DNA_WHERE printf("HERE - %s:%d\n", __FILE__, __LINE__);
#define DNA_WHERE_INT(a)                                \
  printf("HERE - %s:%d = %i\n", __FILE__, __LINE__, a); \
  <
#endif
