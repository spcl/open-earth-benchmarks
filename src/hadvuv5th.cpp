#include <cmath>

// define the domain size and the halo width
int64_t domain_size = 64;
int64_t domain_height = 60;
int64_t halo_width = 4;

typedef double ElementType;

#include "util.h"
#include "hadvuv5th.h"

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

  const int64_t size1D = domain_size + 2 * halo_width;
  const std::array<int64_t, 3> sizes3D = { domain_size + 2*halo_width,
                                           domain_size + 2*halo_width,
                                           domain_height + 2*halo_width };

  Storage3D uin = allocateStorage(sizes3D);
  Storage3D vin = allocateStorage(sizes3D);
  Storage3D uout = allocateStorage(sizes3D);
  Storage3D vout = allocateStorage(sizes3D);
  Storage1D acrlat0 = allocateStorage(size1D);
  Storage1D acrlat1 = allocateStorage(size1D);
  Storage1D tgrlatda0 = allocateStorage(size1D);
  Storage1D tgrlatda1 = allocateStorage(size1D);
  Storage3D uatupos = allocateStorage(sizes3D);
  Storage3D vatupos = allocateStorage(sizes3D);
  Storage3D uatvpos = allocateStorage(sizes3D);
  Storage3D vatvpos = allocateStorage(sizes3D);
  Storage3D uavg = allocateStorage(sizes3D);
  Storage3D vavg = allocateStorage(sizes3D);
  Storage3D ures = allocateStorage(sizes3D);
  Storage3D vres = allocateStorage(sizes3D);


  ElementType eddlat = ldexpl(1.0, -11);
  ElementType eddlon = ldexpl(1.5, -11);

  fillMath(1.0, 3.3, 1.5, 1.5, 2.0, 4.0, uin, domain_size, domain_height);
  fillMath(1.1, 2.0, 1.5, 2.8, 2.0, 4.1, vin, domain_size, domain_height);

  fillMath(1.1, 2.0, 1.5, 2.8, 2.0, 4.1, acrlat0, domain_size, domain_height);
  fillMath(8.0, 9.4, 1.5, 1.7, 2.0, 3.5, acrlat1, domain_size, domain_height);
  fillMath(8.0, 9.4, 1.5, 1.7, 2.0, 3.5, tgrlatda0, domain_size, domain_height);
  fillMath(5.0, 8.0, 1.5, 7.1, 2.0, 4.3, tgrlatda1, domain_size, domain_height);

  initValue(uout, 0.0, domain_size, domain_height);
  initValue(vout, 0.0, domain_size, domain_height);

  hadvuv5th(uout, vout, uin, vin, acrlat0, acrlat1, tgrlatda0, tgrlatda1, uatupos, vatupos, uatvpos, vatvpos, uavg, vavg, ures, vres, eddlat, eddlon);

  // free the storage
  freeStorage(uin);
  freeStorage(vin);
  freeStorage(uout);
  freeStorage(vout);
  freeStorage(acrlat0);
  freeStorage(acrlat1);
  freeStorage(tgrlatda0);
  freeStorage(tgrlatda1);
  freeStorage(uatupos);
  freeStorage(vatupos);
  freeStorage(uatvpos);
  freeStorage(vatvpos);
  freeStorage(uavg);
  freeStorage(vavg);
  freeStorage(ures);
  freeStorage(vres);

  return 0;
}
