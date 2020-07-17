#ifndef HDIFFSMAG_H
#define HDIFFSMAG_H

void hdiffsmag(Storage3D& uout, Storage3D& vout, const Storage3D& uin, const Storage3D& vin, const Storage3D& mask, const Storage1D& crlavo, const Storage1D& crlavu, const Storage1D& crlato, const Storage1D& crlatu, const Storage1D& acrlat0, Storage3D& T_sqr_s, Storage3D& S_sqr_uv, const ElementType eddlat, const ElementType eddlon, const ElementType tau_smag, const ElementType weight_smag) {

  for (int64_t i = 0; i < domain_size; i++) {
    for (int64_t j = 0; j < domain_size; j++) {
      for (int64_t k = 0; k < domain_height; k++) {
        const ElementType frac_1_dx = acrlat0(j) * eddlon;
        const ElementType frac_1_dy = eddlat * EARTH_RADIUS_RECIP;

        // Tension
        const ElementType T_s = (vin(i, j-1, k) - vin(i, j, k)) * frac_1_dy
            - (uin(i-1, j, k) - uin(i, j, k)) * frac_1_dx;
        T_sqr_s(i, j, k) = T_s * T_s; // valid on [1:I,1:J,0:K]
      }
    }
  }

  for (int64_t i = 0; i < domain_size; i++) {
    for (int64_t j = 0; j < domain_size; j++) {
      for (int64_t k = 0; k < domain_height; k++) {
        const ElementType frac_1_dx = acrlat0(j) * eddlon;
        const ElementType frac_1_dy = eddlat * EARTH_RADIUS_RECIP;

        const ElementType S_uv = (uin(i, j+1, k) - uin(i, j, k)) * frac_1_dy
            + (vin(i+1, j, k) - vin(i, j, k)) * frac_1_dx;
        S_sqr_uv(i, j, k) = S_uv * S_uv; // valid on [0:I-1,0:J-1,0:K]
      }
    }
  }

  for (int64_t i = 0; i < domain_size; i++) {
    for (int64_t j = 0; j < domain_size; j++) {
      for (int64_t k = 0; k < domain_height; k++) {
        const ElementType hdweight = weight_smag * mask(i, j, k);

        // I direction
        // valid on [1:I-1,1:J-1,O:K]
        ElementType smag_u = tau_smag * std::sqrt(
            ElementType(0.5) * (T_sqr_s(i+1, j, k) + T_sqr_s(i, j, k))
                + ElementType(0.5) * (S_sqr_uv(i, j-1, k) + S_sqr_uv(i, j, k))
          ) - hdweight;
        smag_u = std::min(ElementType(0.5), std::max(ElementType(0.0), smag_u));

        // valid on [1:I-1,1:J-1,0:K]
        const ElementType lapu = uin(i+1, j, k) + uin(i-1, j, k) - ElementType(2.0) * uin(i, j, k)
            + crlato(j) * (uin(i, j+1, k) - uin(i, j, k))
            + crlatu(j) * (uin(i, j-1, k) - uin(i, j, k));

        // valid on [1:I-1,1:J-1,0:K]
        uout(i, j, k) = uin(i, j, k) + smag_u * lapu;

        // J direction
        // valid on [1:I-1,1:J-1,0:K]
        ElementType smag_v = tau_smag * std::sqrt(
            ElementType(0.5) * (T_sqr_s(i, j+1, k) + T_sqr_s(i, j, k))
                + ElementType(0.5) * (S_sqr_uv(i-1, j, k) + S_sqr_uv(i, j, k))
          ) - hdweight;
        smag_v = std::min(ElementType(0.5), std::max(ElementType(0.0), smag_v));

        // valid on [1:I-1,1:J-1,0:K]
        const ElementType lapv = vin(i+1, j, k) + vin(i-1, j, k) - ElementType(2.0) * vin(i, j, k)
            + crlavo(j) * (vin(i, j+1, k) - vin(i, j, k))
            + crlavu(j) * (vin(i, j-1, k) - vin(i, j, k));

        // valid on [1:I-1,1:J-1,0:K]
        vout(i, j, k) = vin(i, j, k) + smag_v * lapv;
      }
    }
  }
}

#endif // HDIFFSMAG_H
