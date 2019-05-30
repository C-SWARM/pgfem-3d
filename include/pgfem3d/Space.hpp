#pragma once

namespace pgfem3d
{
  /// Defines an enumeration that clients use to query different index spaces.
  enum Space : bool {
    LOCAL  = false,
    GLOBAL = true
  };
}
