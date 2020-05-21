// Common Functions for KMM Spectroscopy Programs - (c) Aleksey Kocherzhenko, 2012

typedef struct {double re; double im;} complex_t;

// complex number operations
complex_t prod(complex_t c1, complex_t c2);
complex_t add(complex_t c1, complex_t c2);
complex_t prod(complex_t c, double d);
complex_t add(complex_t c, double d);
complex_t inv(complex_t c);
complex_t neg(complex_t c);
complex_t conj(complex_t c);
complex_t csqrt(complex_t c); 
complex_t cexp(complex_t c);

// Calculate THETA(k) = exp[-2*i*theta(k)]
complex_t getTH(double k, double COUP);

// Calculate allowed k values for the "+" branch (alpha) and "-" branch (beta)
int AlphaBetaVals(double *pkplus, double *pkminus, int UNITS);

// Calculate energy of single-exciton state with given k
double getEn(double k, double COUP);

// Function to print the binary representation of a number
void printbin(unsigned int n);

// Function to check whether the number n has a 1 at bit position m
unsigned int checkbin(unsigned int n, unsigned int m);

