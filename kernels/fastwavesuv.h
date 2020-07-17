#ifndef FASTWAVESUV_H
#define FASTWAVESUV_H

void fastwavesuv(Storage3D& uout, Storage3D& vout, const Storage3D& uin, const Storage3D& vin, const Storage3D& utens, const Storage3D& vtens, const Storage3D& wgtfac, const Storage3D& ppuv, const Storage3D& hhl, const Storage3D& rho, const Storage1D& fx, Storage3D& ppgk, Storage3D& ppgc, Storage3D& ppgu, Storage3D& ppgv, const ElementType edadlat, const ElementType dt) {

  for (int64_t i = 0; i < domain_size + 1; ++i) {
    for (int64_t j = 0; j < domain_size + 1; ++j) {
      for (int64_t k = 0; k < domain_height + 1; ++k) {
        ppgk(i, j, k) = wgtfac(i, j, k) * ppuv(i, j, k) +
                        (ElementType(1.0) - wgtfac(i, j, k)) * ppuv(i, j, k - 1);
      }
    }
  }

  for (int64_t i = 0; i < domain_size + 1; ++i) {
    for (int64_t j = 0; j < domain_size + 1; ++j) {
      for (int64_t k = 0; k < domain_height; ++k) {
        ppgc(i, j, k) = ppgk(i, j, k + 1) - ppgk(i, j, k);
      }
    }
  }

  for (int64_t i = 0; i < domain_size; ++i) {
    for (int64_t j = 0; j < domain_size; ++j) {
      for (int64_t k = 0; k < domain_height; ++k) {
        ppgu(i, j, k) = (ppuv(i + 1, j, k) - ppuv(i, j, k)) +
                        (ppgc(i + 1, j, k) + ppgc(i, j, k)) * ElementType(0.5) *
                            ((hhl(i, j, k + 1) + hhl(i, j, k)) -
                             (hhl(i + 1, j, k + 1) + hhl(i + 1, j, k))) /
                            ((hhl(i, j, k + 1) - hhl(i, j, k)) +
                             (hhl(i + 1, j, k + 1) - hhl(i + 1, j, k)));
      }
    }
  }

  for (int64_t i = 0; i < domain_size; ++i) {
    for (int64_t j = 0; j < domain_size; ++j) {
      for (int64_t k = 0; k < domain_height; ++k) {
        ppgv(i, j, k) = (ppuv(i, j + 1, k) - ppuv(i, j, k)) +
                        (ppgc(i, j + 1, k) + ppgc(i, j, k)) * ElementType(0.5) *
                            ((hhl(i, j, k + 1) + hhl(i, j, k)) -
                             (hhl(i, j + 1, k + 1) + hhl(i, j + 1, k))) /
                            ((hhl(i, j, k + 1) - hhl(i, j, k)) +
                             (hhl(i, j + 1, k + 1) - hhl(i, j + 1, k)));
      }
    }
  }

  for (int64_t i = 0; i < domain_size; ++i) {
    for (int64_t j = 0; j < domain_size; ++j) {
      for (int64_t k = 0; k < domain_height; ++k) {
        uout(i, j, k) =
            uin(i, j, k) +
            dt * (utens(i, j, k) -
                    ppgu(i, j, k) * (ElementType(2.0)*fx(j) / (rho(i + 1, j, k) + rho(i, j, k))));
      }
    }
  }

  for (int64_t i = 0; i < domain_size; ++i) {
    for (int64_t j = 0; j < domain_size; ++j) {
      for (int64_t k = 0; k < domain_height; ++k) {
        vout(i, j, k) =
            vin(i, j, k) +
            dt * (vtens(i, j, k) -
                    ppgv(i, j, k) * (ElementType(2.0)*edadlat / (rho(i, j + 1, k) + rho(i, j, k))));

      }
    }
  }
}

#endif // FASTWAVESUV_H
