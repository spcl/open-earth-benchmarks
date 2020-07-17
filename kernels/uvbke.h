#ifndef UVBKE_H
#define UVBKE_H

void uvbke(Storage3D& ub, Storage3D& vb, const Storage3D& uc, const Storage3D& vc, const Storage3D& cosa, const Storage3D& rsina) {
  for (int64_t k = 0; k < domain_height; ++k) {
    for (int64_t i = 0; i < domain_size; ++i) {
      for (int64_t j = 0; j < domain_size; ++j) {
        ub(i, j, k) =
            ((dt5 *
              ((uc(i, j - 1, k) + uc(i, j, k)) -
               ((vc(i - 1, j, k) + vc(i, j, k)) *
                cosa(i, j, k)))) * rsina(i, j, k));
      }
    }
  }

  for (int64_t k = 0; k < domain_height; ++k) {
    for (int64_t i = 0; i < domain_size; ++i) {
      for (int64_t j = 0; j < domain_size; ++j) {
        vb(i, j, k) =
            ((dt5 *
              ((vc(i - 1, j, k) + vc(i, j, k)) -
               ((uc(i, j - 1, k) + uc(i, j, k)) *
                cosa(i, j, k)))) * rsina(i, j, k));
      }
    }
  }
}

#endif // UVBKE_H
