#ifndef NH_P_GRAD_H
#define NH_P_GRAD_H

void nh_p_grad(Storage3D& uout, Storage3D& vout, const Storage3D& uin, const Storage3D& vin, const Storage3D& rdx, const Storage3D& rdy, const Storage3D& gz, const Storage3D& pp, const Storage3D& pk3, const Storage3D& wk1, Storage3D& wk, Storage3D& du, Storage3D& dv, const ElementType dt) {

  for (int64_t k = 0; k < domain_height; ++k) {
    for (int64_t i = 0; i < domain_size+1; ++i) {
      for (int64_t j = 0; j < domain_size+1; ++j) {
        wk(i, j, k) = (pk3(i, j, k + 1) - pk3(i, j, k));
      }
    }
  }

  for (int64_t k = 0; k < domain_height; ++k) {
    for (int64_t i = 0; i < domain_size; ++i) {
      for (int64_t j = 0; j < domain_size; ++j) {
        du(i, j, k) =
            ((dt / (wk(i, j, k) + wk(i + 1, j, k))) *
             (((gz(i, j, k + 1) - gz(i + 1, j, k)) * (pk3(i + 1, j, k + 1) - pk3(i, j, k))) +
              ((gz(i, j, k) - gz(i + 1, j, k + 1)) * (pk3(i, j, k + 1) - pk3(i + 1, j, k)))));
        uout(i, j, k) =
            (((uin(i, j, k) + du(i, j, k)) +
              ((dt / (wk1(i, j, k) + wk1(i + 1, j, k))) *
               (((gz(i, j, k + 1) - gz(i + 1, j, k)) * (pp(i + 1, j, k + 1) - pp(i, j, k))) +
                ((gz(i, j, k) - gz(i + 1, j, k + 1)) * (pp(i, j, k + 1) - pp(i + 1, j, k)))))) *
             rdx(i, j, k));
      }
    }
  }

  for (int64_t k = 0; k < domain_height; ++k) {
    for (int64_t i = 0; i < domain_size; ++i) {
      for (int64_t j = 0; j < domain_size; ++j) {
        dv(i, j, k) =
            ((dt / (wk(i, j, k) + wk(i, j + 1, k))) *
             (((gz(i, j, k + 1) - gz(i, j + 1, k)) * (pk3(i, j + 1, k + 1) - pk3(i, j, k))) +
              ((gz(i, j, k) - gz(i, j + 1, k + 1)) * (pk3(i, j, k + 1) - pk3(i, j + 1, k)))));
        vout(i, j, k) =
            (((vin(i, j, k) + dv(i, j, k)) +
              ((dt / (wk1(i, j, k) + wk1(i, j + 1, k))) *
               (((gz(i, j, k + 1) - gz(i, j + 1, k)) * (pp(i, j + 1, k + 1) - pp(i, j, k))) +
                ((gz(i, j, k) - gz(i, j + 1, k + 1)) * (pp(i, j, k + 1) - pp(i, j + 1, k)))))) *
             rdy(i, j, k));
      }
    }
  }
}

#endif // NH_P_GRAD_H
