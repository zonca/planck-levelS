#include "levels_facilities.h"
#include "error_handling.h"

int main (int argc, const char **argv)
  {
  PLANCK_DIAGNOSIS_BEGIN
  alm2grid_module (argc, argv);
  PLANCK_DIAGNOSIS_END
  }
