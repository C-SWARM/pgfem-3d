#pragma once

#include <algorithm>
#include <cassert>
#include <type_traits>
#include <vector>

namespace pgfem3d {
template <class T>
class JaggedArray {
 public:
  /// Allocate an empty jagged array.
  JaggedArray() : data_(), offsets_() {
  }

  /// Allocate a jagged array with row widths given in cols.
  ///
  /// @tparam         U A numeric type that stores the row widths.
  ///
  /// @param       cols An array with the column width for each row.
  template <class U>
  JaggedArray(const std::vector<U>& cols) : JaggedArray() {
    static_assert(std::is_arithmetic<U>::value, "cols must be a numeric array");
    resize(cols);
  }

  /// Release all of the memory associated with the Array.
  void release() {
    offsets_.resize(0);
    offsets_.shrink_to_fit();
    data_.resize(0);
    data_.shrink_to_fit();
  }

  /// Resize the array to have `rows` row where the `op` provides the columns in
  /// each row.
  ///
  /// @tparam        Op An operator type that provides the column count of each
  ///                   row.
  ///
  /// @param       rows The number of rows to resize the array to.
  /// @param         op The operator that provides the number of columns in each
  ///                   row.
  template <class Op>
  void resize(size_t rows, Op&& op) {
    if (rows == 0) {
      release();
      return;
    }

    offsets_.resize(rows);
    offsets_.shrink_to_fit();
    offsets_[0] = 0;
    for (size_t i = 1, e = rows; i < e; ++i) {
      offsets_[i] = offsets_[i - 1] + op(i - 1);
    }
    data_.resize(offsets_[rows - 1] + op(rows - 1));
    data_.shrink_to_fit();
  }

  /// Resize the array using the array `cols` as the new shape.
  ///
  /// @tparam         U A numeric type that stores the row widths.
  ///
  /// @param       cols An array of row widths.
  template <class U>
  void resize(const std::vector<U>& cols) {
    static_assert(std::is_arithmetic<U>::value, "cols must be a numeric array");
    resize(cols.size(), [&](auto i) {
      return cols[i];
    });
  }

  /// Treat array as a 2d array.
  ///
  /// @param          i The row index.
  /// @param          j The column index.
  ///
  /// @returns          The value at `i,j`.
  const T& operator()(size_t i, size_t j) const {
    assert(i < rows() and j < cols(i));
    return data_[offsets_[i] + j];
  }

  /// Treat array as a 2d array.
  ///
  /// @param          i The row index.
  /// @param          j The column index.
  ///
  /// @returns          A reference to the value at `i,j`.
  T& operator()(size_t i, size_t j) {
    assert(i < rows() and j < cols(i));
    return data_[offsets_[i] + j];
  }

  template <class U>
  T& append(size_t i, std::vector<U>& next) {
    assert(next.size() == rows());
    return this->operator()(i, next[i]++);
  }

  /// Get the number of rows in the array.
  size_t rows() const {
    return offsets_.size();
  }

  /// Get the number of columns in a row of the array.
  size_t cols(size_t i) const {
    assert(i < rows());
    if (i == rows() - 1) {
      return data_.size() - offsets_[i];
    }
    else {
      return offsets_[i + 1] - offsets_[i];
    }
  }

  /// Sort and compact the array to eliminate duplicates in rows.
  void compact() {
    size_t insert = 0;
    for (size_t i = 0, e = rows(); i < e; ++i) {
      sort(i);
      offsets_[i] = compact(i, insert);
    }
    data_.resize(insert);
    data_.shrink_to_fit();
  }

 private:
  void sort(size_t i) {
    assert(cols(i) > 0);
    std::sort(&data_[offsets_[i]], &data_[offsets_[i]] + cols(i));
  }

  /// Compact a sorted row, using `insert` as the location to insert into.
  ///
  /// This will update the insert position as the compaction happens. It will
  /// return the original insertion location which is the offset to the
  /// compacted row.
  size_t compact(size_t i, size_t& insert) {
    auto at = [this](auto i, auto j) -> T& {
      return this->operator()(i, j);
    };

    auto offset = insert;
    data_[insert++] = at(i, 0);
    for (size_t j = 1, e = cols(i); j < e; ++j) {
      if (at(i, j - 1) < at(i, j)) {
        data_[insert++] = at(i, j);
      }
    }
    return offset;
  }

  std::vector<T> data_;                  //!< The array data stored contiguously
  std::vector<size_t> offsets_;          //!< The row offsets.
}; // class JaggedArray<T>
} // namespace pgfem3d
