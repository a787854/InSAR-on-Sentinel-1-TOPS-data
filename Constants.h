#ifndef CONSTANTS_H
#define  CONSTANTS_H

#include <math.h>

static char eachLine[4 * 258];//For reading a line
#define SOL 299792458.0
#define PI 3.1415926535897931e+0
#define Interval 2047
#define max(a,b)  ((a)>(b)?(a):(b))
#define min(a,b)  ((a)<(b)?(a):(b))
#define demNoData -32768.0

inline double sqr(const double &x) { return(x*x); }
inline double sinc(const double &x){ return ((x == 0) ? 1 : sin(PI*x) / (PI*x)); }

inline double max4(double a, double b, double c, double d)
{
	return max(max(a, b), max(c, d));
}

inline double min4(double a, double b, double c, double d)
{
	return min(min(a, b), min(c, d));
}

inline double rect(const double &x)
{
	double ans = 0.0;
	if (x<0.5 && x>-0.5) ans = 1;
	else if (x == 0.5 || x == -0.5) ans = 0.5;
	return ans;
}

inline double deg2rad(const double x)
{
	return x * PI / 180.0;
}

inline double abs3D(double x[3]){ return sqrt(sqr(x[0]) + sqr(x[1]) + sqr(x[2])); }



#endif//CONSTANTS_H