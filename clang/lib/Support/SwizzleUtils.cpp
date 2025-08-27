#include <algorithm>
#include <clang/Support/SwizzleUtils.h>
#include <optional>

namespace clang {
namespace swizzle {

std::unordered_map<swizzle_key, std::optional<swizzle_solver::solution_type>>
    swizzle_solver::cache;

swizzle_key::swizzle_key(int _bits, const std::vector<access_pattern> &_pattern)
    : bits{_bits}, hash{}, patterns{} {
  for (auto i : _pattern) {
    std::sort(i.begin(), i.end());
    patterns.push_back(i);
    hash ^= std::hash<access_pattern>{}(i);
  }
  std::sort(patterns.begin(), patterns.end(),
            [](const access_pattern &l, const access_pattern &r) {
              return std::hash<access_pattern>{}(l) <
                     std::hash<access_pattern>{}(r);
            });
}

bool swizzle_key::operator==(const swizzle_key &_other) const noexcept {
  if (_other.hash != hash)
    return false;
  if (_other.patterns.size() != patterns.size())
    return false;
  for (auto i = 0u; i < patterns.size(); i++) {
    if (patterns[i] != _other.patterns[i])
      return false;
  }
  return true;
}

std::optional<swizzle_solver::solution_type>
swizzle_solver::solve(int bits,
                      const std::vector<access_pattern> &patterns) noexcept {
  const auto bit_mask = (1u << bits) - 1;
  
  for(const auto& i: patterns){
    if(!check_pattern(i)) return std::optional<solution_type>{};
  }

  std::vector<ker_space_basis> w; // kernal space vectors.
  // {v, d} means vector keeps dim 'd' as primary dimension and use it to
  // eliminate dimension of any other vectors

  // maintaining patterns to be elimiated by vectors from w
  std::vector<access_pattern> elim_patterns{patterns};

  // record used primary dimensions
  auto prohibited_bits = 0u;

  // dim w should be at least bits - 3 to construct a 3 dimension complement
  // space
  for (auto collected_dims = 0; collected_dims < bits - 3; collected_dims++) {
    //const auto remain_bits = bits - collected_dims;

    // allowed dimensions to be primary dimensions for next vector
    const auto allow_bits = bit_mask & ~prohibited_bits;

    auto step_solved = false;

    // iterate over subsets of allow_bits. ivec != allow_bits
    for (auto ivec = (allow_bits - 1) & allow_bits;;
         ivec = (ivec - 1) & allow_bits) {
      // next vector to be added into w
      const auto vec = allow_bits & ~ivec;
      // selected primary dimension of vec
      const auto dim_bit = __builtin_ffs(vec) - 1;

      // try to use vec to eliminate patterns.
      // if patterns keep their ranks, the w \cap {vec} cannot represent any
      // vector from the image space of P_i
      auto rank_lost = false;
      std::vector<access_pattern> next_elim_patterns{elim_patterns};
      for (auto &p : next_elim_patterns) {
        for (auto i = 0; i < 3; i++) {
          if ((p[i] >> dim_bit) & 0x1)
            p[i] ^= vec;
        }

        if (!check_pattern(p)) {
          rank_lost = true;
          break;
        }
      }
      if (!rank_lost) { // vec is acceptable
        elim_patterns = std::move(next_elim_patterns);
        prohibited_bits |= 1u << dim_bit;
        w.push_back(std::make_tuple(vec, dim_bit));
        step_solved = true;
        break;
      }
      if (!ivec)
        break;
    }
    // required step cannot be resolved so no solution
    if (!step_solved)
      return std::optional<solution_type>{};
  }

  // simplify the w space to make sure all non-free dimension depend only on
  // free dimensions.
  for (auto i = w.rbegin(); i != w.rend(); i++) {
    auto &current_vec = *i;
    for (auto j = i + 1; j != w.rend(); j++) {
      if (std::get<0>(*j) & (1u << std::get<1>(current_vec)))
        *j = std::make_tuple(std::get<0>(*j) ^ std::get<0>(current_vec),
                             std::get<1>(*j));
    }
  }

  std::vector<uint32_t> basic_solution;
  uint32_t free_dims = bit_mask & ~prohibited_bits;
  while (free_dims) {
    auto lowbit = free_dims & (-free_dims);
    auto result = lowbit;

    free_dims ^= lowbit;

    for (auto &[vec, dim] : w) {
      if (vec & lowbit)
        result |= 1u << dim;
    }
    // generate complement basis vector as row vectors of basic solution
    basic_solution.push_back(result);
  }

  int min_diag = bits + 3;
  std::array<uint32_t, 3> opt_ans;

  // iterate over GL(3,2) to generate better solution
  for (auto i = 1; i < 8; i++) {
    for (auto j = i + 1; j < 8; j++) {
      for (auto k = j + 1; k < 8; k++) {
        if ((i ^ j) == k)
          continue;

        std::array<int, 3> perm{i, j, k};

        do {
          std::array<uint32_t, 3> candidate{gemv(basic_solution, perm[0]),
                                            gemv(basic_solution, perm[1]),
                                            gemv(basic_solution, perm[2])};
          auto diags = count_diag(bits, candidate);
          if (min_diag > diags) {
            opt_ans = candidate;
            min_diag = diags;
          }
        } while (std::next_permutation(perm.begin(), perm.end()));
      }
    }
  }

  return std::optional<solution_type>{{opt_ans, min_diag}};
}

} // namespace swizzle
} // namespace clang