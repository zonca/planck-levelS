#ifndef ZOD_INTERFACE_H
#define ZOD_INTERFACE_H

#ifdef __cplusplus
extern "C" {
#endif

void zle_init(double freq, const int *component, const double *fact, int ncomp,
  double accuracy, int nloop, double mindist, double maxdist);

void zle_compute(const double *sun, const double *bary, int nsamples,
  const double *directions, double *signal);

void zle_destroy(void);

#ifdef __cplusplus
}
#endif

#endif
