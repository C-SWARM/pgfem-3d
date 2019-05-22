#pragma once

#include <algorithm>                            // std::upper_bound
#include <cassert>                              // assert
#include <numeric>                              // std::partial_sum
#include <vector>                               // std::vector

namespace pgfem3d {
/// This template represents a blocked row distribution.
///
/// Blocked row distributions are single dimensional distributions of a
/// contiguous set of rows across ranks. They support operations like checking
/// to see if a rank owns a row, which rank owns a row, and the minimum and
/// maximum row for a rank.
///
/// They are used in the PGFem3D communication infrastructure to distribute the
/// global sparse system for assembly.
template <class RowId, class RankId = int>
class BlockedRowDistribution {
 public:
  /// Allocate an empty distribution.
  BlockedRowDistribution() : offsets_(1) {
    offsets_[0] = 0;
  }

  /// Construct a blocked row distribution for a number of ranks.
  ///
  /// This constructor requires the number of ranks and an array of rows per
  /// rank.
  ///
  ///
  BlockedRowDistribution(RankId nRanks, const RowId blockSizes[])
      : BlockedRowDistribution()
  {
    rebuild(nRanks, blockSizes);
  }

  /// Rebuild a blocked row distribution for a number of ranks.
  void rebuild(RankId nRanks, const RowId blockSizes[]) {
    offsets_.resize(nRanks + 1);
    offsets_[0] = 0;
    std::partial_sum(&blockSizes[0], &blockSizes[nRanks], &offsets_[1]);
  }

  /// Lookup the owner of a row.
  RankId owner(RowId i) const {
    assert(0 <= i and i < rows());
    auto begin = std::begin(offsets_);
    auto   end = std::end(offsets_);
    auto    it = std::upper_bound(begin, end, i);
    auto     n = std::distance(begin, it);
    assert(0 < n and n < ranks() + 1);
    return RankId(n - 1);
  }

  /// Check if a row is owned by a rank.
  bool owns(RankId rank, RowId i) const {
    assert(0 <= i and i < rows());
    return (min(rank) <= i and i < max(rank));
  }

  /// Get the rank's minimum row id.
  RowId min(RankId rank) const {
    assert(0 <= rank and rank < ranks());
    return offsets_[rank];
  }

  /// Get the rank's maximum row id.
  RowId max(RankId rank) const {
    assert(0 <= rank and rank < ranks());
    return offsets_[rank + 1];
  }

  /// Get the number of ranks in the distribution.
  RankId ranks() const {
    size_t n = offsets_.size();
    assert(0 < n and n < size_t(std::numeric_limits<RankId>::max()));
    return RankId(n) - 1;
  }

  /// Get the number of rows in the distribution.
  RowId rows() const {
    return offsets_.back();
  }

  /// Get the number of rows for a particular rank.
  auto rows(RankId rank) const {
    return max(rank) - min(rank);
  }

 private:
  std::vector<RowId> offsets_;
};
}
