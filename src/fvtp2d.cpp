#include <cmath>

// define the domain size and the halo width
int64_t domain_size = 64;
int64_t domain_height = 60;
int64_t halo_width = 4;

typedef double ElementType;

#define p1 ElementType(0.583333)  // 7/12
#define p2 ElementType(0.083333)  // 1/12

#include "util.h"
#include "fvtp2d.h"

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
  Storage3D q    = allocateStorage(sizes3D);
  Storage3D crx  = allocateStorage(sizes3D);
  Storage3D cry  = allocateStorage(sizes3D);
  Storage3D ra_x = allocateStorage(sizes3D);
  Storage3D ra_y = allocateStorage(sizes3D);
  Storage3D xfx  = allocateStorage(sizes3D);
  Storage3D yfx  = allocateStorage(sizes3D);
  Storage3D area = allocateStorage(sizes3D);

  Storage3D q_i = allocateStorage(sizes3D);
  Storage3D q_j = allocateStorage(sizes3D);
  Storage3D fx1 = allocateStorage(sizes3D);
  Storage3D fx2 = allocateStorage(sizes3D);
  Storage3D fy1 = allocateStorage(sizes3D);
  Storage3D fy2 = allocateStorage(sizes3D);

  Storage3D fxx  = allocateStorage(sizes3D);
  Storage3D fyy  = allocateStorage(sizes3D);
  Storage3D al = allocateStorage(sizes3D);
  Storage3D almq = allocateStorage(sizes3D);
  Storage3D br   = allocateStorage(sizes3D);
  Storage3D b0   = allocateStorage(sizes3D);
  Storage3D smt5 = allocateStorage(sizes3D);

  fillMath(1.0, 3.3, 1.5, 1.5, 2.0, 4.0, q, domain_size, domain_height);
  fillMath(1.1, 2.0, 1.5, 2.8, 2.0, 4.1, crx, domain_size, domain_height);
  fillMath(1.1, 2.0, 1.5, 2.8, 2.0, 4.1, cry, domain_size, domain_height);
  fillMath(8.0, 9.4, 1.5, 1.7, 2.0, 3.5, ra_x, domain_size, domain_height);
  fillMath(8.0, 9.4, 1.5, 1.7, 2.0, 3.5, ra_y, domain_size, domain_height);
  fillMath(5.0, 8.0, 1.5, 7.1, 2.0, 4.3, xfx, domain_size, domain_height);
  fillMath(5.0, 8.0, 1.5, 7.1, 2.0, 4.3, yfx, domain_size, domain_height);
  fillMath(3.2, 7.0, 2.5, 6.1, 3.0, 2.3, area, domain_size, domain_height);

  initValue(q_i, -1.0, domain_size, domain_height);
  initValue(q_j, -1.0, domain_size, domain_height);
  initValue(fx1, -1.0, domain_size, domain_height);
  initValue(fx2, -1.0, domain_size, domain_height);
  initValue(fy1, -1.0, domain_size, domain_height);
  initValue(fy2, -1.0, domain_size, domain_height);

  fvtp2d(q_i, q_j, fx1, fx2, fy1, fy2, q, crx, cry, ra_x, ra_y, xfx, yfx, area, fxx, fyy, al, almq, br, b0, smt5);

  // free the storage
  freeStorage(q);
  freeStorage(crx);
  freeStorage(cry);
  freeStorage(ra_x);
  freeStorage(ra_y);
  freeStorage(xfx);
  freeStorage(yfx);
  freeStorage(area);

  freeStorage(q_i);
  freeStorage(q_j);
  freeStorage(fx1);
  freeStorage(fx2);
  freeStorage(fy1);
  freeStorage(fy2);

  freeStorage(fxx);
  freeStorage(fyy);
  freeStorage(al);
  freeStorage(almq);
  freeStorage(br);
  freeStorage(b0);
  freeStorage(smt5);

  return 0;
}
