double r8_abs ( double x );
double r8_epsilon ( );
void r8_fehl ( void f ( double t, double y[], double yp[], double params[] ), int neqn,
  double y[], double t, double h, double yp[], double f1[], double f2[], double f3[], 
  double f4[], double f5[], double s[] , double params[]);
double r8_max ( double x, double y );
double r8_min ( double x, double y );
int r8_rkf45 ( void f ( double t, double y[], double yp[], double params[] ), int neqn,
  double y[], double yp[], double *t, double tout, double *relerr, double abserr, 
  int flag, double params[] );
double r8_sign ( double x );

void timestamp ( );
