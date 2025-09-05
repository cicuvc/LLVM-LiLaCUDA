#ifndef _CLANG_SUPPORT_SWIZZLEUTILS_H
#define _CLANG_SUPPORT_SWIZZLEUTILS_H

#include <array>
#include <cstdint>
#include <optional>
#include <unordered_map>
#include <vector>

namespace clang {
namespace swizzle {

using AccessPattern = std::array<uint32_t, 3>;
using KerSpaceBasis = std::tuple<uint32_t, int>;

struct SwizzleSolver;
struct SwizzleKey;

struct SwizzleKey {
  int bits;
  int hash;
  std::vector<AccessPattern> patterns;

  SwizzleKey(int _bits, const std::vector<AccessPattern> &_pattern);

  bool operator==(const SwizzleKey &_other) const noexcept;
};

} // namespace swizzle
} // namespace clang

template <> struct std::hash<clang::swizzle::SwizzleKey> {
  inline std::size_t
  operator()(clang::swizzle::SwizzleKey const &s) const noexcept {
    return s.hash;
  }
};

template <> struct std::hash<clang::swizzle::AccessPattern> {
  std::size_t
  operator()(clang::swizzle::AccessPattern const &s) const noexcept {
    return s[0] ^ s[1] ^ s[2];
  }
};

namespace clang {
namespace swizzle {

struct SwizzleSolver {
  using SolutionType = std::tuple<std::array<uint32_t, 3>, int>;

  static std::unordered_map<SwizzleKey, std::optional<SolutionType>> cache;

  // Check if rank of access pattern == 3
  static inline bool checkPattern(const AccessPattern &ap) noexcept {
    return (ap[0] != ap[1]) && (ap[0] != ap[2]) && (ap[1] != ap[2]) &&
           ((ap[0] ^ ap[1]) != ap[2]) && (ap[0]) && ap[1] && ap[2];
  }

  // Count non-zero diagonal of swizzle matrix
  static inline int countDiag(int bits,
                               const std::array<uint32_t, 3> &swizzle) noexcept {

    auto result = 1;
    result += !!(swizzle[2] & 0x1);
    result += ((swizzle[1] & 0x1)) || ((swizzle[2] >> 1) & 0x1);
    for (auto i = 1; i < bits; i++) {
      result += ((swizzle[0] >> i) & 0x1) || ((swizzle[1] >> (i + 1)) & 0x1) ||
                ((swizzle[2] >> (i + 2)) & 0x1);
    }
    return result;
  }

  // calculate F_2 matrix left gemv
  static inline uint32_t gemv(const std::vector<uint32_t> &results, int vec) noexcept {
    auto result = 0u;
    for (auto i = results.begin(); i != results.end(); i++, vec >>= 1) {
      if (vec & 0x1)
        result ^= *i;
    }
    return result;
  }

  // main procedure to solve swizzle matrix for given patterns
  static std::optional<SolutionType>
  solve(int bits, const std::vector<AccessPattern> &patterns) noexcept;

  static inline std::optional<SolutionType>
  getSwizzleSolution(int _bits,
                       const std::vector<AccessPattern> &_pattern) noexcept {
    SwizzleKey key{_bits, _pattern};
    const auto it = cache.find(key);
    if (it != cache.end())
      return it->second;
    const auto result = solve(_bits, _pattern);
    cache.emplace(key, result);
    return result;
  }
};
} // namespace swizzle
} // namespace clang

#endif