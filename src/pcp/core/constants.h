#pragma once

namespace pcp {

// Portable math constants (M_PI is not available on MSVC without extra macros).
inline constexpr double kPi = 3.14159265358979323846;
inline constexpr double kDegToRad = kPi / 180.0;
inline constexpr double kRadToDeg = 180.0 / kPi;

}  // namespace pcp
