#ifndef HDIFFSA_H
#define HDIFFSA_H

void hdiffsa(Storage3D& out, const Storage3D& in, const Storage3D& mask, const Storage1D& crlato, const Storage1D& crlatu, Storage3D& lap, Storage3D& flx, Storage3D& fly) {

  for(int64_t k=0; k < domain_height; ++k) {
    for(int64_t i=-1; i < domain_size+1; ++i) {
      for(int64_t j=-1; j < domain_size+1; ++j) {
        lap(i,j,k) = in(i-1,j,k) + in(i+1,j,k) - ElementType(2.0)*in(i,j,k) +
                crlato(j) * (in(i,j+1,k) - in(i,j,k)) +
                crlatu(j) * (in(i,j-1,k) - in(i,j,k));
      }
    }
    for(int64_t i=-1; i < domain_size+1; ++i) {
      for(int64_t j=-1; j < domain_size+1; ++j) {
        flx(i,j,k) = lap(i+1,j,k) - lap(i,j,k);

        if(flx(i,j,k) * (in(i+1,j,k) - in(i,j,k)) > 0) {
          flx(i,j,k) = 0;
        }

        fly(i,j,k) = crlato(j)*(lap(i,j+1,k) - lap(i,j,k));

        if(fly(i,j,k) * (in(i,j+1,k) - in(i,j,k)) > 0) {
          fly(i,j,k) = 0;
        }
      }
    }
    for(int64_t i=0; i < domain_size; ++i) {
      for(int64_t j=0; j < domain_size; ++j) {
        out(i,j,k) = in(i,j,k) + (flx(i-1,j,k) - flx(i,j,k) +
                fly(i,j-1,k) - fly(i,j,k))*mask(i,j,k);
      }
    }
  }
}

#endif // HDIFFSA_H
