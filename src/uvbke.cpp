#include <cmath>

// define the domain size and the halo width
int64_t domain_size = 64;
int64_t domain_height = 60;
int64_t halo_width = 4;

typedef double ElementType;

const double dt = 225.0;
const double dt5 = 0.5 * dt;

#include "util.h"
#include "uvbke.h"

// program times the execution of the linked program and times the result
int main(int argc, char **argv) {

  if(argc == 3) {
      domain_size = atoi(argv[1]);
      domain_height = atoi(argv[2]);
  } else if (argc == 1) {
  } else {
      std::cout << "Either provide the domain size and domain height like this \"./kernel 128 60\" or do not provide any arguments at all in which case the program is ran with domain size 64 and domain heigh 60" << std::endl;
      exit(1);
  }

  const std::array<int64_t, 3> sizes3D = { domain_size + 2*halo_width,
                                           domain_size + 2*halo_width,
                                           domain_height + 2*halo_width };


  // allocate the storage
  Storage3D uc  = allocateStorage(sizes3D);
  Storage3D vc  = allocateStorage(sizes3D);
  Storage3D cosa  = allocateStorage(sizes3D);
  Storage3D rsina  = allocateStorage(sizes3D);

  Storage3D ub = allocateStorage(sizes3D);
  Storage3D vb = allocateStorage(sizes3D);

  fillMath(1.0, 3.3, 1.5, 1.5, 2.0, 4.0, uc, domain_size, domain_height);
  fillMath(1.1, 2.0, 1.5, 2.8, 2.0, 4.1, vc, domain_size, domain_height);
  fillMath(4.0, 1.7, 1.5, 6.3, 2.0, 1.4, cosa, domain_size, domain_height);
  fillMath(8.0, 9.4, 1.5, 1.7, 2.0, 3.5, rsina, domain_size, domain_height);

  initValue(ub, -1.0, domain_size, domain_height);
  initValue(vb, -1.0, domain_size, domain_height);

  uvbke(ub, vb, uc, vc, cosa, rsina);

  // free the storage
  freeStorage(uc);
  freeStorage(vc);
  freeStorage(cosa);
  freeStorage(rsina);

  freeStorage(ub);
  freeStorage(vb);

  return 0;
}
