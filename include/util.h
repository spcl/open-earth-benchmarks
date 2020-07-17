#ifndef UTIL_H
#define UTIL_H

#include <iostream>
#include <array>

#define EARTH_RADIUS ((ElementType)6371.229e3) // radius of the earth
#define EARTH_RADIUS_RECIP ((ElementType)1.0 / EARTH_RADIUS)

const ElementType pi(std::acos(-1.0));

namespace {

template <typename ElementType, size_t Rank>
struct Storage {
    ElementType *allocatedPtr;
    ElementType *alignedPtr;
    int64_t offset;
    std::array<int64_t, Rank> sizes;   // omitted when rank == 0
    std::array<int64_t, Rank> strides; // omitted when rank == 0

   template<typename... T>
   struct type { };

   template <typename... T>
   ElementType &operator()(T... arg) {
       return operator()(type<T...>(), arg...);
   }

   template <typename... T>
   const ElementType &operator()(T... arg) const {
       return operator()(type<T...>(), arg...);
   }

   template<typename... T>
   ElementType &operator()(type<T...>, T... arg);

   template<typename... T>
   const ElementType &operator()(type<T...>, T... arg) const;
};

typedef Storage<ElementType, 1> Storage1D;
typedef Storage<ElementType, 2> Storage2D;
typedef Storage<ElementType, 3> Storage3D;

template <> template <>
ElementType &Storage<ElementType, 1>::operator()<int64_t>(Storage::type<int64_t>, int64_t i) {
  return alignedPtr[offset + strides[0] * i];
}

template <> template <>
const ElementType &Storage<ElementType, 1>::operator()<int64_t>(Storage::type<int64_t>, int64_t i) const {
  return alignedPtr[offset + strides[0] * i];
}

template <> template <>
ElementType &Storage<ElementType, 2>::operator()<int64_t, int64_t>(Storage::type<int64_t, int64_t>, int64_t i, int64_t j) {
    return alignedPtr[offset + strides[1] * i + strides[0] * j];
}

template <> template <>
const ElementType &Storage<ElementType, 2>::operator()<int64_t, int64_t>(Storage::type<int64_t, int64_t>, int64_t i, int64_t j) const {
    return alignedPtr[offset + strides[1] * i + strides[0] * j];
}

template <> template <>
ElementType &Storage<ElementType, 3>::operator()<int64_t, int64_t, int64_t>(Storage::type<int64_t, int64_t, int64_t>, int64_t i, int64_t j, int64_t k) {
    return alignedPtr[offset + strides[2] * i + strides[1] * j + strides[0] * k];
}

template <> template <>
const ElementType &Storage<ElementType, 3>::operator()<int64_t, int64_t, int64_t>(Storage::type<int64_t, int64_t, int64_t>, int64_t i, int64_t j, int64_t k) const {
    return alignedPtr[offset + strides[2] * i + strides[1] * j + strides[0] * k];
}

} // namespace

// allocate Storages
Storage1D allocateStorage(const int64_t size) {
  Storage1D result;
  // initialize the size
  result.sizes[0] = size;
  // initialize the strides
  result.strides[0] = 1;
  result.offset = halo_width * result.strides[0];
  result.allocatedPtr = new ElementType[size + (32 - halo_width)];
  result.alignedPtr = &result.allocatedPtr[(32 - halo_width)];
  return result;
}

Storage2D allocateStorage(const std::array<int64_t, 2> sizes) {
  Storage2D result;
  // initialize the size
  result.sizes[1] = sizes[0];
  result.sizes[0] = sizes[1];
  // initialize the strides
  result.strides[1] = 1;
  result.strides[0] = result.sizes[1];
  result.offset = halo_width * result.strides[0] +
                  halo_width * result.strides[1];
  const int64_t allocSize = sizes[0] * sizes[1];
  result.allocatedPtr = new ElementType[allocSize + (32 - halo_width)];
  result.alignedPtr = &result.allocatedPtr[(32 - halo_width)];
  return result;
}

Storage3D allocateStorage(const std::array<int64_t, 3> sizes) {
  Storage3D result;
  // initialize the size
  result.sizes[2] = sizes[0];
  result.sizes[1] = sizes[1];
  result.sizes[0] = sizes[2];
  // initialize the strides
  result.strides[2] = 1;
  result.strides[1] = result.sizes[2];
  result.strides[0] = result.sizes[2] * result.sizes[1];
  result.offset = halo_width * result.strides[0] +
                  halo_width * result.strides[1] +
                  halo_width * result.strides[2];
  const int64_t allocSize = sizes[0] * sizes[1] * sizes[2];
  result.allocatedPtr = new ElementType[allocSize + (32 - halo_width)];
  result.alignedPtr = &result.allocatedPtr[(32 - halo_width)];
  return result;
}

template <typename Storage>
void freeStorage(Storage &ref) {
  delete ref.allocatedPtr;
  ref.allocatedPtr = nullptr;
  ref.alignedPtr = nullptr;
}

void fillMath(ElementType a, ElementType b, ElementType c, ElementType d,
              ElementType e, ElementType f, Storage3D &field, const int64_t domain_size, const int64_t domain_height) {
  ElementType dx = ElementType(1.0) / (ElementType) (domain_size + 2*halo_width);
  ElementType dy = ElementType(1.0) / (ElementType) (domain_size + 2*halo_width);

  for (int64_t j = -halo_width; j < domain_size + halo_width; j++) {
    for (int64_t i = -halo_width; i < domain_size + halo_width; i++) {
      ElementType x = dx * (ElementType) i;
      ElementType y = dy * (ElementType) j;
      for (int64_t k = 0; k < domain_height; k++) {
          field(i, j, k) = k*ElementType(10e-3) + a*(b + cos(pi*(x + c*y)) + sin(d*pi*(x + e*y)))/f;
      }
    }
  }
}

void fillMath(ElementType a, ElementType b, ElementType c, ElementType d,
              ElementType e, ElementType f, Storage2D &field, const int64_t domain_size, const int64_t domain_heigh) {
  ElementType dx = ElementType(1.0) / (ElementType) (domain_size + 2*halo_width);
  ElementType dy = ElementType(1.0) / (ElementType) (domain_size + 2*halo_width);

  for (int64_t j = -halo_width; j < domain_size + halo_width; j++) {
    for (int64_t i = -halo_width; i < domain_size + halo_width; i++) {
      ElementType x = dx * (ElementType) i;
      ElementType y = dy * (ElementType) j;
      field(i, j) = a*(b + cos(pi*(x + c*y)) + sin(d*pi*(x + e*y)))/f;
    }
  }
}

void fillMath(ElementType a, ElementType b, ElementType c, ElementType d,
              ElementType e, ElementType f, Storage1D &field, const int64_t domain_size, const int64_t domain_heigh) {
  ElementType dx = ElementType(1.0) / (ElementType) (domain_size + 2*halo_width);

  for (int64_t i = -halo_width; i < domain_size + halo_width; i++) {
    ElementType x = dx * (ElementType) i;
    field(i) = a*(b + cos(pi*(c*x)) + sin(d*pi*(e*x)))/f;
  }
}

void initValue(Storage3D &ref, const ElementType val, const int64_t domain_size, const int64_t domain_heigh) {
  for (int64_t i = -halo_width; i < domain_size + halo_width; ++i)
    for (int64_t j = -halo_width; j < domain_size + halo_width; ++j)
      for (int64_t k = -halo_width; k < domain_height + halo_width; ++k) {
        ref(i, j, k) = val;
      }
}

#endif // UTIL_H
