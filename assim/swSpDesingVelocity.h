#pragma once
/**
 * @file swSpDesingVelocity.h
 * @brief Kurganov-Petrova desingularized velocity, shared by the native and the
 *        differentiated (_E) MUSCL-EVR reconstructions.
 *
 * Added for automatic differentiation; not part of native SW2D.
 *
 * This lives in a header rather than being duplicated in swSpHypMusclEvr.cpp and
 * swSpHypMusclEvr_E.cpp on purpose: both reconstructions MUST use the exact same
 * regularization. A divergence between the classic and the _E path here would
 * silently linearize the adjoint around a different trajectory.
 */

#include <cmath>
#include <algorithm>

namespace swSpDesing {

// Kurganov-Petrova desingularized velocity u = q / h.
// For h >= h_desing it returns EXACTLY q/h (no change to the physics in wet cells);
// for h < h_desing it transitions smoothly to 0 with a BOUNDED derivative.
// This matters for the Enzyme adjoint: the raw q/max(h,small_real) has
// d(q/h)/dh = -q/h^2 ~ -1e9 near dry cells, which makes the reverse-mode
// adjoint state ||d_h|| explode over the assimilation window. The desingularized
// form keeps the wave speeds and reconstructed discharges differentiable.
//   h_desing : depth scale below which the velocity is regularized (tunable).
inline double desingVelocity(double q, double h, double h_desing)
{
    double h2 = h * h;
    double h4 = h2 * h2;
    double e4 = h_desing * h_desing * h_desing * h_desing;
    double denom = std::sqrt(h4 + std::max(h4, e4));
    return 1.4142135623730951 * h * q / denom;   // sqrt(2) * h * q / sqrt(h^4 + max(h^4, e^4))
}

// Regularization depth scale (m). >= a few mm bounds d(u)/dh enough to tame the
// adjoint explosion while leaving cells deeper than h_desing physically untouched.
constexpr double kHDesing = 5e-3;

} // namespace swSpDesing
