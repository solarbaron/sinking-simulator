// SPDX-License-Identifier: MIT
//
// The optics of a two-zone fire, with no device in sight.
//
// Everything `shaders/smoke.frag` computes is stated in `engine/gpu/smoke.cpp` as
// well, and this file checks *that* statement against closed forms:
// Stefan-Boltzmann for the whole-spectrum integral, Wien for where its peak is,
// Beer-Lambert for the transmittance, and plain geometry for the two path
// lengths. `tests/test_smoke_render.cpp` then holds the GPU against the same
// functions, so a drift between the shader and the model shows up as a pixel
// disagreement rather than as two implementations quietly agreeing on the wrong
// thing.
#include "engine/gpu/smoke.hpp"
#include "engine/sim/fire.hpp"
#include "game/prototype/ferry.hpp"
#include "harness.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

using testing::expectEqual;
using testing::expectNear;
using testing::expectTrue;

namespace {

// --- Planck ---------------------------------------------------------------

// The whole spectrum, integrated, must be `sigma T^4 / pi`. That is the one
// closed form the quadrature has no way to fake: it fixes the constant, the
// exponent and the normalisation at once, and it is sensitive to the whole
// range rather than to the visible band the renderer actually uses.
void testTheSpectrumIntegratesToStefanBoltzmann() {
    // Derived from h, c and k rather than quoted -- the same discipline
    // `fire.hpp` applies to the caloric constants. CODATA's 5.670374419e-8 is
    // this expression evaluated at the 2019 defining constants and then **rounded
    // to ten significant digits**, so the two agree to ten and differ in the
    // eleventh: 3.2525e-11 relative, measured. That gap is the published rounding
    // rather than an error here, so 3.3e-11 is the tightest this may honestly be.
    //
    // The comment here used to say "exactly", and read the agreement as a check on
    // the algebra at 1e-12 -- a tolerance the code never used and which the real
    // gap fails by 33x. What the tightening does buy is rejecting a *nine*-digit
    // quotation: 5.67037442e-8 is 1.44e-10 off, which 1e-9 accepted.
    expectNear("Stefan-Boltzmann comes out of h, c and k", gpu::stefanBoltzmann(),
               5.670374419e-8, 5.670374419e-8 * 3.3e-11);

    // The same constant under the other convention: `fire.hpp` quotes CODATA's
    // ten-digit figure, this file derives it from the SI's exact h, c and k. Two
    // computations of one physical constant, and nothing in the suite related them
    // -- `test_fire.cpp` holds the quoted one against its own local re-derivation,
    // so a divergence between the *modules* would have shown up nowhere and left
    // the fire model and the renderer radiating at quietly different rates.
    // Allowed to differ, but only by that published rounding.
    expectNear("and it is fire.hpp's kStefanBoltzmann under the other convention",
               gpu::stefanBoltzmann(), sim::fire::kStefanBoltzmann,
               sim::fire::kStefanBoltzmann * 3.3e-11);

    // 1 nm to 1 m. **The upper limit is not decoration.** Stopping at 1 mm loses
    // 5.6e-6 of a 300 K body's power to the Rayleigh-Jeans tail -- measured, by
    // running exactly this check with that limit -- which is eight orders of
    // magnitude above the tolerance below and would have been read as a broken
    // quadrature rather than as a truncated integral.
    //
    // **5e-14 and not 1e-9.** The worst of these four residuals is 6.7e-15, and
    // moving from 512 to 2048 intervals takes it to 2.5e-14, so 5e-14 bounds the
    // rule's own rounding floor rather than one lucky configuration. It matters
    // beyond tidiness: at 1e-9 the whole suite stayed green with `stefanBoltzmann()`
    // replaced by the quoted literal -- measured, by substituting it -- so the
    // derivation this file exists to justify was unfalsifiable. At 5e-14 all four
    // of these go red against the literal, which is what the header claims.
    for (double temperature : {300.0, 800.0, 1200.0, 2000.0}) {
        const double integrated = gpu::blackbodyBandRadiance(1e-9, 1.0, temperature, 512);
        const double closed = gpu::blackbodyRadiance(temperature);
        expectNear("the integrated spectrum is sigma T^4 / pi at " +
                       std::to_string(static_cast<int>(temperature)) + " K",
                   integrated / closed, 1.0, 5e-14);
    }

    // And the quadrature converges rather than happening to land: refining it must
    // reduce the error, and the coarse rule must be *visibly* worse. Without this
    // second half the check above passes on a rule that is exact by accident at
    // one node count.
    const double closed = gpu::blackbodyRadiance(1200.0);
    const double coarse = std::abs(gpu::blackbodyBandRadiance(1e-9, 1.0, 1200.0, 4) / closed - 1.0);
    const double fine = std::abs(gpu::blackbodyBandRadiance(1e-9, 1.0, 1200.0, 512) / closed - 1.0);
    expectTrue("the coarse rule is measurably wrong", coarse > 1e-4);
    expectTrue("and refining it fixes that", fine < 1e-6 * coarse);
}

// Planck's law must reduce to Rayleigh-Jeans, `B -> 2 c k T / lambda^4`, as
// `x = hc / (lambda k T)` goes to zero -- and the *rate* is the test, not a value
// at one wavelength. The leading correction is `x / 2`, so the agreement has to
// improve linearly as x shrinks.
//
// That is what catches a `exp(x) - 1` written where `expm1(x)` belongs: the
// cancellation there costs a relative error of about eps/x, which *grows* as the
// limit is approached, so the same sweep gets worse instead of better. It is worth
// stating because the obvious check does not catch it at all -- the whole-spectrum
// integral above moves by 1e-16 between the two forms, since the region that
// cancels carries almost none of the power.
void testPlanckReachesItsRayleighJeansLimit() {
    const double temperature = 300.0;
    for (double x : {1e-6, 1e-8, 1e-10}) {
        const double wavelength = gpu::kPlanck * gpu::kLightSpeed /
                                  (x * gpu::kBoltzmann * temperature);
        const double rayleighJeans = 2.0 * gpu::kLightSpeed * gpu::kBoltzmann * temperature /
                                     (wavelength * wavelength * wavelength * wavelength);
        const double relative =
            std::abs(gpu::planckRadiance(wavelength, temperature) / rayleighJeans - 1.0);
        // The correction is x/2; asserting against x leaves a factor of two of
        // room and no more, which is what makes this a rate rather than a bound.
        //
        // `std::to_string` on a double is six decimal places, which prints every
        // x here as 0.000000 and every failure as the same message. Scientific
        // notation, so the line names the case it failed on.
        char what[160];
        std::snprintf(what, sizeof what,
                      "Planck reaches Rayleigh-Jeans as hc/lambda kT falls to %.0e:"
                      " relative %.3e", x, relative);
        expectTrue(what, relative < x);
    }
    // Vacuity: the limit has to be somewhere the formula is *not* already exact,
    // or a constant function would pass the sweep.
    const double visible = gpu::kPlanck * gpu::kLightSpeed / (10.0 * gpu::kBoltzmann * temperature);
    const double rj = 2.0 * gpu::kLightSpeed * gpu::kBoltzmann * temperature /
                      (visible * visible * visible * visible);
    expectTrue("and it is nowhere near it in the band a fire is seen in",
               gpu::planckRadiance(visible, temperature) < 0.01 * rj);
}

// Wien's displacement law, against the published constant rather than against a
// root this file solved for itself.
void testThePeakSitsWhereWienSaysItDoes() {
    constexpr double kWien = 2.897771955e-3;  // m K, CODATA
    for (double temperature : {600.0, 1100.0, 1800.0, 5772.0}) {
        const double predicted = kWien / temperature;
        // Ternary search over two decades around the prediction. Wide enough that
        // finding the prediction is a result rather than an assumption.
        double lo = predicted * 0.1, hi = predicted * 10.0;
        for (int i = 0; i < 200; ++i) {
            const double a = lo + (hi - lo) / 3.0, b = hi - (hi - lo) / 3.0;
            if (gpu::planckRadiance(a, temperature) < gpu::planckRadiance(b, temperature))
                lo = a;
            else
                hi = b;
        }
        expectNear("Planck peaks at Wien's wavelength at " +
                       std::to_string(static_cast<int>(temperature)) + " K",
                   0.5 * (lo + hi) / predicted, 1.0, 1e-6);
    }
}

// The three display bands, and the one display mapping in the file.
void testTheBandsAreASpectralSampleAndTheMappingIsOneNumber() {
    double hot[3], cool[3];
    gpu::blackbodyBands(1400.0, hot);
    gpu::blackbodyBands(700.0, cool);

    expectTrue("the bands are red, green, blue in that order and R is brightest below 4000 K",
               hot[0] > hot[1] && hot[1] > hot[2]);
    expectTrue("every band grows with temperature",
               hot[0] > cool[0] && hot[1] > cool[1] && hot[2] > cool[2]);
    // The visible bands are a sliver of a 1400 K body's output, which is why the
    // upper layer of a fire glows dimly rather than blazing.
    const double visible = hot[0] + hot[1] + hot[2];
    expectTrue("the visible bands are a small part of the whole",
               visible < 0.02 * gpu::blackbodyRadiance(1400.0));
    // Each band is exactly the integral over its own edges: the split is the
    // spectrum's, not a colour transform's.
    expectNear("the red band is the integral over 600-700 nm", hot[0],
               gpu::blackbodyBandRadiance(600e-9, 700e-9, 1400.0, 64), hot[0] * 1e-12);

    // The display rule's interval count is settled here rather than assumed.
    // Against a 4096-interval reference the worst disagreement is at the cold end,
    // where the blue band sits furthest down the Wien tail -- which is why the
    // sweep starts at 300 K and covers all three bands rather than the red one a
    // fire is actually seen in.
    double worst = 0.0;
    for (double temperature : {300.0, 600.0, 900.0, 1400.0, 2500.0})
        for (int band = 0; band < 3; ++band) {
            const double reference = gpu::blackbodyBandRadiance(
                gpu::kBandEdges[band], gpu::kBandEdges[band + 1], temperature, 4096);
            const double cheap = gpu::blackbodyBandRadiance(
                gpu::kBandEdges[band], gpu::kBandEdges[band + 1], temperature,
                gpu::kDisplayBandIntervals);
            worst = std::max(worst, std::abs(cheap / reference - 1.0));
        }
    expectTrue("the display band rule is the 4096-interval answer", worst < 1e-12);
    // And a coarser rule is *not*, so the count is a measurement rather than a
    // round number that happened to work. Two intervals slip to 1e-8 in the blue
    // band at 300 K; one slips at 300 K in the red band too.
    const double coarse = gpu::blackbodyBandRadiance(gpu::kBandEdges[0], gpu::kBandEdges[1],
                                                     300.0, 2);
    const double exact = gpu::blackbodyBandRadiance(gpu::kBandEdges[0], gpu::kBandEdges[1],
                                                    300.0, 4096);
    expectTrue("where a coarser one would not be", std::abs(coarse / exact - 1.0) > 1e-10);

    gpu::SmokeShading shading;
    double reference[3], below[3];
    gpu::emissiveColour(shading.referenceTemperature, shading, reference);
    expectNear("the reference temperature saturates its brightest channel",
               std::max({reference[0], reference[1], reference[2]}), 1.0, 1e-12);
    gpu::emissiveColour(shading.referenceTemperature * 0.8, shading, below);
    expectTrue("and everything cooler is dimmer in every channel",
               below[0] < reference[0] && below[1] < reference[1] && below[2] < reference[2]);

    // The finding worth having a number for: the exponential is brutal, so no
    // linear display mapping shows a 900 K layer and a 1200 K one at once.
    double warm[3], hotter[3];
    gpu::emissiveColour(900.0, shading, warm);
    gpu::emissiveColour(1200.0, shading, hotter);
    expectTrue("900 K to 1200 K is more than two decades of visible radiance",
               hotter[0] / warm[0] > 100.0);

    // Ambient air does not glow. Not a threshold -- Planck does it.
    double ambient[3];
    gpu::emissiveColour(sim::kTAmbient, shading, ambient);
    expectTrue("ambient air emits nothing the display can show", ambient[0] < 1e-12);
}

// --- Beer-Lambert ---------------------------------------------------------

void testTransmittanceIsBeerLambert() {
    expectNear("a slab of k = 0.25 over 4 m transmits exp(-1)",
               gpu::transmittance(0.25, 4.0), std::exp(-1.0), 1e-15);

    // Doubling the path squares the transmittance. In double this is an identity
    // rather than a tolerance, and it is the property that catches an integration
    // that has picked up an extra term.
    for (double k : {0.01, 0.3, 1.7, 12.0}) {
        for (double d : {0.2, 1.0, 3.5}) {
            const double single = gpu::transmittance(k, d);
            const double doubled = gpu::transmittance(k, 2.0 * d);
            expectNear("doubling the path squares the transmittance", doubled, single * single,
                       std::max(1e-16, single * single * 1e-14));
        }
    }

    // Exactly 1.0, not nearly: this is what makes an unsmoked frame bit-identical.
    expectTrue("zero extinction transmits exactly one",
               gpu::transmittance(0.0, 1e6) == 1.0 && gpu::transmittance(1e6, 0.0) == 1.0);
}

// --- The geometry ---------------------------------------------------------

gpu::SmokeVolume unitBox() {
    gpu::SmokeVolume v;
    v.name = "box";
    v.lo = {-2.0, -3.0, 0.0};
    v.hi = {2.0, 3.0, 6.0};
    v.interfaceZ = 4.0;
    return v;
}

void testTheRayMeetsThePrismWhereGeometrySaysItDoes() {
    const gpu::SmokeVolume v = unitBox();
    double tEnter = 0, tExit = 0;

    expectTrue("a ray along -x through the middle hits",
               gpu::intersectVolume(v, {10.0, 0.0, 3.0}, {-1.0, 0.0, 0.0}, 1e30, tEnter, tExit));
    expectNear("and enters at the +x face", tEnter, 8.0, 1e-12);
    expectNear("and leaves at the -x face", tExit, 12.0, 1e-12);

    expectTrue("a ray that passes above it misses",
               !gpu::intersectVolume(v, {10.0, 0.0, 9.0}, {-1.0, 0.0, 0.0}, 1e30, tEnter, tExit));
    expectTrue("a ray pointing away from it misses",
               !gpu::intersectVolume(v, {10.0, 0.0, 3.0}, {1.0, 0.0, 0.0}, 1e30, tEnter, tExit));

    // Started inside: the near limit is the eye, not the box.
    expectTrue("a ray starting inside it hits",
               gpu::intersectVolume(v, {0.0, 0.0, 3.0}, {-1.0, 0.0, 0.0}, 1e30, tEnter, tExit));
    expectNear("from zero", tEnter, 0.0, 1e-15);
    expectNear("to the near face", tExit, 2.0, 1e-12);

    // An opaque surface inside the box truncates the segment rather than
    // rejecting it -- which is the whole reason the depth buffer is sampled and
    // not tested against.
    expectTrue("an opaque surface inside the box truncates the segment",
               gpu::intersectVolume(v, {10.0, 0.0, 3.0}, {-1.0, 0.0, 0.0}, 9.5, tEnter, tExit));
    expectNear("at the surface", tExit, 9.5, 1e-12);
    expectTrue("and an opaque surface in front of it hides it entirely",
               !gpu::intersectVolume(v, {10.0, 0.0, 3.0}, {-1.0, 0.0, 0.0}, 5.0, tEnter, tExit));

    // Exactly parallel to a face and inside it: the substituted reciprocal has to
    // give the same answer an infinity would, without ever meeting a zero.
    expectTrue("a ray parallel to a face still hits",
               gpu::intersectVolume(v, {0.0, 10.0, 3.0}, {0.0, -1.0, 0.0}, 1e30, tEnter, tExit));
    expectNear("entering at the +y face", tEnter, 7.0, 1e-12);
    expectNear("and leaving at the -y face", tExit, 13.0, 1e-12);
}

void testThePathSplitsAtTheInterface() {
    const gpu::SmokeVolume v = unitBox();  // floor 0, deckhead 6, interface at 4
    double tEnter = 0, tExit = 0;

    // Straight down the middle from above: 2 m of upper layer, then 4 m of lower.
    expectTrue("a vertical ray hits",
               gpu::intersectVolume(v, {0.0, 0.0, 20.0}, {0.0, 0.0, -1.0}, 1e30, tEnter, tExit));
    gpu::LayerPath path = gpu::layerPath(v, {0.0, 0.0, 20.0}, {0.0, 0.0, -1.0}, tEnter, tExit);
    expectNear("2 m of hot layer", path.upper, 2.0, 1e-12);
    expectNear("4 m of cool layer", path.lower, 4.0, 1e-12);
    expectTrue("and the descending ray meets the hot layer first", path.upperFirst);

    // The same ray from below: same lengths, opposite order. If the order were
    // taken from where the eye is rather than from the ray's slope, this would
    // agree with the case above and be wrong.
    expectTrue("a ray from below hits",
               gpu::intersectVolume(v, {0.0, 0.0, -20.0}, {0.0, 0.0, 1.0}, 1e30, tEnter, tExit));
    path = gpu::layerPath(v, {0.0, 0.0, -20.0}, {0.0, 0.0, 1.0}, tEnter, tExit);
    expectNear("still 2 m of hot layer", path.upper, 2.0, 1e-12);
    expectNear("still 4 m of cool layer", path.lower, 4.0, 1e-12);
    expectTrue("but the cool layer comes first", !path.upperFirst);

    // Horizontal, wholly inside one layer.
    expectTrue("a horizontal ray in the hot layer hits",
               gpu::intersectVolume(v, {10.0, 0.0, 5.0}, {-1.0, 0.0, 0.0}, 1e30, tEnter, tExit));
    path = gpu::layerPath(v, {10.0, 0.0, 5.0}, {-1.0, 0.0, 0.0}, tEnter, tExit);
    expectNear("all of it is hot", path.upper, 4.0, 1e-12);
    expectNear("and none of it cool", path.lower, 0.0, 1e-15);

    expectTrue("a horizontal ray in the cool layer hits",
               gpu::intersectVolume(v, {10.0, 0.0, 1.0}, {-1.0, 0.0, 0.0}, 1e30, tEnter, tExit));
    path = gpu::layerPath(v, {10.0, 0.0, 1.0}, {-1.0, 0.0, 0.0}, tEnter, tExit);
    expectNear("all of it is cool", path.lower, 4.0, 1e-12);
    expectNear("and none of it hot", path.upper, 0.0, 1e-15);

    // A slanted ray, where the split is neither zero nor everything. 45 degrees
    // down through the interface: the crossing is at a known height and the two
    // lengths are the geometry, not a fraction of the span.
    const sim::Vec3 down = normalize(sim::Vec3{-1.0, 0.0, -1.0});
    expectTrue("a slanted ray hits",
               gpu::intersectVolume(v, {10.0, 0.0, 11.0}, down, 1e30, tEnter, tExit));
    path = gpu::layerPath(v, {10.0, 0.0, 11.0}, down, tEnter, tExit);
    // Entry at x = 2 after 8 m of x travel, so z = 11 - 8 = 3: already below the
    // interface. It then leaves through the *floor* at z = 0, which is x = -1 and
    // not the far wall at x = -2 -- so the segment is three metres of x travel,
    // not four. The first version of this expectation said four; the code was
    // right and the arithmetic in the test was wrong, which is the order CLAUDE.md
    // asks the question in.
    expectNear("the slanted ray is wholly in the cool layer", path.lower,
               3.0 * std::sqrt(2.0), 1e-12);
    expectNear("with nothing hot", path.upper, 0.0, 1e-15);

    // And one that genuinely straddles: enter at z = 5, cross at z = 4 one metre
    // of drop later.
    expectTrue("a straddling ray hits",
               gpu::intersectVolume(v, {10.0, 0.0, 13.0}, down, 1e30, tEnter, tExit));
    path = gpu::layerPath(v, {10.0, 0.0, 13.0}, down, tEnter, tExit);
    expectNear("one metre of drop in the hot layer", path.upper, std::sqrt(2.0), 1e-12);
    expectNear("and three in the cool one", path.lower, 3.0 * std::sqrt(2.0), 1e-12);
    expectTrue("hot first", path.upperFirst);
    expectNear("and the two account for the whole segment", path.upper + path.lower,
               tExit - tEnter, 1e-12);
}

// A heeled ship. `fire.hpp` says the interface is horizontal in the *body* frame,
// so the whole medium has to rotate with her -- and the path lengths through it
// must be exactly what they were, because a rotation is an isometry. Getting the
// transform's handedness wrong produces a picture that is plausible from every
// single viewpoint.
void testTheMediumRotatesWithTheShip() {
    gpu::SmokeVolume upright = unitBox();
    const sim::Vec3 eye{10.0, 0.0, 13.0};
    const sim::Vec3 direction = normalize(sim::Vec3{-1.0, 0.0, -1.0});
    double t0 = 0, t1 = 0;
    expectTrue("the upright case hits",
               gpu::intersectVolume(upright, eye, direction, 1e30, t0, t1));
    const gpu::LayerPath reference = gpu::layerPath(upright, eye, direction, t0, t1);

    // 25 degrees of heel about the bow axis, with the ray carried round with her.
    const sim::Mat3 heel = sim::Quat::fromAxisAngle({1, 0, 0}, 25.0 * sim::kDegToRad).toMat3();
    gpu::SmokeVolume heeled = upright;
    heeled.rotation = heel;
    heeled.translation = {3.0, -1.0, 0.5};
    const sim::Vec3 movedEye = heel * eye + heeled.translation;
    const sim::Vec3 movedDirection = heel * direction;

    double h0 = 0, h1 = 0;
    expectTrue("the heeled case hits too",
               gpu::intersectVolume(heeled, movedEye, movedDirection, 1e30, h0, h1));
    const gpu::LayerPath moved = gpu::layerPath(heeled, movedEye, movedDirection, h0, h1);
    expectNear("a rotation does not change the hot path", moved.upper, reference.upper, 1e-12);
    expectNear("nor the cool one", moved.lower, reference.lower, 1e-12);
    expectTrue("nor which is met first", moved.upperFirst == reference.upperFirst);
    // Vacuity: the reference must have had both layers in it, or a renderer that
    // lost the interface entirely would pass this.
    expectTrue("and the reference had both layers in it",
               reference.upper > 0.5 && reference.lower > 0.5);
}

// --- The transfer integral -------------------------------------------------

void testTheCompositeIsTheTransferEquation() {
    gpu::SmokeVolume v = unitBox();
    const double background[3] = {0.2, 0.4, 0.8};
    double out[3];

    // A clear volume leaves the background alone, bit for bit. This is the CPU
    // half of the zero-smoke identity the GPU test asserts on pixels.
    gpu::LayerPath path;
    path.upper = 3.0;
    path.lower = 2.0;
    v.upper.emission[0] = v.upper.emission[1] = v.upper.emission[2] = 5.0;
    gpu::compositeOver(v, path, background, out);
    expectTrue("a clear but glowing-hot volume changes nothing at all",
               out[0] == background[0] && out[1] == background[1] && out[2] == background[2]);

    // One absorbing layer over a background: Beer-Lambert exactly.
    v.upper = {};
    v.upper.extinction = 0.5;
    path.lower = 0.0;
    gpu::compositeOver(v, path, background, out);
    for (int c = 0; c < 3; ++c)
        expectNear("one absorbing layer is exp(-k d) times the background", out[c],
                   background[c] * std::exp(-1.5), 1e-15);

    // Two layers, hot one in front: written out by hand rather than by calling the
    // same expression a second time.
    v.upper.extinction = 0.5;
    v.upper.emission[0] = 0.9;
    v.lower.extinction = 0.2;
    v.lower.emission[0] = 0.1;
    path.upper = 3.0;
    path.lower = 2.0;
    path.upperFirst = true;
    gpu::compositeOver(v, path, background, out);
    const double tu = std::exp(-0.5 * 3.0), tl = std::exp(-0.2 * 2.0);
    const double byHand = 0.9 * (1.0 - tu) + tu * (0.1 * (1.0 - tl) + tl * background[0]);
    expectNear("the two-slab composite is the transfer equation", out[0], byHand, 1e-15);

    // Order matters, and it matters because of emission. Swapping which layer the
    // ray meets first must change the answer -- and must stop changing it when
    // neither layer glows.
    double swapped[3];
    path.upperFirst = false;
    gpu::compositeOver(v, path, background, swapped);
    expectTrue("which layer is met first changes the answer",
               std::abs(swapped[0] - out[0]) > 1e-3);
    expectNear("but not the transmittance", swapped[2], out[2], 1e-15);

    v.upper.emission[0] = 0.0;
    v.lower.emission[0] = 0.0;
    double darkFirst[3], darkSecond[3];
    path.upperFirst = true;
    gpu::compositeOver(v, path, background, darkFirst);
    path.upperFirst = false;
    gpu::compositeOver(v, path, background, darkSecond);
    expectNear("two dark layers compose in either order", darkSecond[0], darkFirst[0], 1e-15);
}

// --- The camera -------------------------------------------------------------

void testTheDepthBasisComesOutOfTheMatrix() {
    const double nearPlane = 0.7, farPlane = 900.0;
    const sim::Vec3 eye{40.0, -25.0, 18.0};
    const sim::Vec3 target{2.0, 3.0, 5.0};
    const sim::Mat4 mvp = sim::perspective(48.0 * sim::kDegToRad, 16.0 / 9.0, nearPlane, farPlane) *
                          sim::lookAt(eye, target, {0, 0, 1});
    float matrix[16];
    mvp.toFloats(matrix);

    gpu::DepthBasis basis;
    expectTrue("a perspective matrix yields a depth basis", gpu::depthBasisFrom(matrix, basis));
    const sim::Vec3 forward = normalize(target - eye);
    expectNear("the view direction comes back out of it", dot(basis.forward, forward), 1.0, 1e-7);
    expectNear("and so does a = far / (near - far)", basis.a, farPlane / (nearPlane - farPlane),
               1e-6);
    expectNear("and b = near far / (near - far)", basis.b,
               nearPlane * farPlane / (nearPlane - farPlane), 1e-3);

    // The round trip, which is what the shader actually does: project a known
    // point, take its depth, and recover the distance along the view axis.
    for (const sim::Vec3& probe :
         {sim::Vec3{0, 0, 0}, sim::Vec3{-10, 12, 3}, sim::Vec3{5, -30, 25}}) {
        double clip[4];
        mvp.transform(probe, clip);
        const double depth = clip[2] / clip[3];
        const double recovered = basis.b / (depth + basis.a);
        expectNear("depth inverts to the distance along the view axis", recovered,
                   dot(forward, probe - eye), 1e-6 * std::abs(dot(forward, probe - eye)));
    }

    // The clear value is the far plane, which is what makes an empty pixel stop
    // the ray at the frustum rather than at zero.
    //
    // **This is the worst-conditioned point of the inversion and the tolerance says
    // so.** `a` is about -1.0008 for these planes, so `1 + a` is a cancellation
    // down to 7.8e-4, and the matrix arrives as float32: seven digits in `a` leave
    // four in `1 + a`. Measured at 8e-6 relative, asserted at 1e-4. It costs
    // nothing, because the far plane is where nothing is -- every ray that reaches
    // it has already left the medium.
    expectNear("a cleared depth of 1 is the far plane", basis.b / (1.0 + basis.a), farPlane,
               farPlane * 1e-4);

    // An orthographic projection has no centre of projection, and the renderer
    // refuses it rather than reconstructing a distance that does not exist.
    float orthographic[16] = {1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1};
    expectTrue("an orthographic matrix is refused",
               !gpu::depthBasisFrom(orthographic, basis));
}

// --- The ship ---------------------------------------------------------------

sim::fire::Model ferryFire(sim::Ship& ship, const sim::Sea& sea, int steps) {
    sim::fire::Model model;
    model.attach(ship, {ship.findCompartment("engine_room_s"),
                        ship.findCompartment("engine_room_p")});
    sim::fire::DesignFire d;
    d.name = "machinery";
    d.compartment = model.findGas("engine_room_s");
    d.baseZ = 2.5;
    d.diameter = 2.5;
    d.growthCoefficient = sim::fire::kGrowthFast;
    d.peakHeatRelease = 4.0e6;
    d.steadyDuration = 600.0;
    model.fires.push_back(d);
    for (int i = 0; i < steps; ++i) model.step(1.0, ship, sea);
    return model;
}

// Both layers, on a compartment built by hand with tracer in each of them.
//
// **The ferry cannot ask this question.** Its lower layer never carries any
// products at all -- the plume draws *from* the lower layer and deposits into the
// upper one, and on this ship nothing puts tracer back down -- so every check on a
// real fire compares a lower extinction of exactly zero against another exactly
// zero, and a renderer that read the upper layer's volume for both would pass all
// of them. That is a hole a mutation pass finds and a scenario does not.
void testBothLayersTakeTheirOwnConcentration() {
    sim::fire::GasCompartment g;
    g.floorZ = 0.0;
    g.ceilingZ = 10.0;
    g.floorArea = 40.0;
    g.gasVolume = 400.0;
    // Deliberately different in every term, so a swapped index cannot cancel.
    g.upper.mass = 30.0;
    g.upper.energy = 30.0 * sim::fire::kCvAir * 900.0;
    g.upper.products = 0.9;
    g.lower.mass = 300.0;
    g.lower.energy = 300.0 * sim::fire::kCvAir * 320.0;
    g.lower.products = 0.06;

    const double upperVolume = g.upperVolume();
    const double lowerVolume = g.gasVolume - upperVolume;
    expectTrue("the fixture splits the space unevenly",
               upperVolume > 1.0 && lowerVolume > 1.0 && upperVolume < 0.5 * lowerVolume);

    gpu::SmokeShading shading;
    const gpu::SmokeLayer upper = gpu::layerFrom(g.upper, upperVolume, shading);
    const gpu::SmokeLayer lower = gpu::layerFrom(g.lower, lowerVolume, shading);
    expectNear("the hot layer's extinction is its own products over its own volume",
               upper.extinction, shading.massExtinction * g.upper.products / upperVolume,
               upper.extinction * 1e-12);
    expectNear("and the cool layer's is its own", lower.extinction,
               shading.massExtinction * g.lower.products / lowerVolume,
               lower.extinction * 1e-12);
    // The two differ by more than any rounding, so reading one volume for both is
    // a failure rather than a coincidence.
    expectTrue("and the two are far apart", upper.extinction > 20.0 * lower.extinction);

    double hot[3], cool[3];
    gpu::emissiveColour(900.0, shading, hot);
    gpu::emissiveColour(320.0, shading, cool);
    expectNear("each layer emits Planck at its own temperature", upper.emission[0], hot[0],
               hot[0] * 1e-12);
    expectNear("including the cool one", lower.emission[0], cool[0],
               std::max(1e-30, cool[0] * 1e-12));
    expectTrue("which at 320 K is nothing at all", cool[0] < 1e-20 && hot[0] > 1e-3);

    // A layer with no soot in it does not glow, however hot: Kirchhoff, and the
    // reason a clear volume is bit-identical to no volume at all.
    sim::fire::Layer clear = g.upper;
    clear.products = 0.0;
    const gpu::SmokeLayer none = gpu::layerFrom(clear, upperVolume, shading);
    expectTrue("a sootless layer has no extinction", none.extinction == 0.0);
    expectNear("though it is still at 900 K", none.emission[0], hot[0], hot[0] * 1e-12);
}

// The prism drawn is the prism the model solved on. `Model::attach` says in as
// many words that a bounding box would over-state the floor area "by the turn of
// the bilge"; drawing the bounding box would put that error back.
void testTheDrawnVolumeIsTheModelsVolume() {
    sim::Ship ship = game::buildFerry();
    const sim::Sea sea;
    ship.initialise(sea);
    const sim::fire::Model model = ferryFire(ship, sea, 300);
    const gpu::SmokeShading shading;
    const std::vector<gpu::SmokeVolume> volumes = gpu::volumesFromFire(model, ship, shading);

    expectEqual("every tracked gas space is drawable", static_cast<long long>(volumes.size()),
                static_cast<long long>(model.gas.size()));

    for (std::size_t i = 0; i < volumes.size(); ++i) {
        const gpu::SmokeVolume& v = volumes[i];
        const sim::fire::GasCompartment& g = model.gas[i];
        expectNear("the drawn plan area is the model's floor area", v.planArea(), g.floorArea,
                   g.floorArea * 1e-12);
        const double drawn = v.planArea() * (v.hi.z - v.lo.z);
        expectNear("so the drawn gas volume is the model's gas volume", drawn, g.gasVolume,
                   g.gasVolume * 1e-12);
        const double drawnUpper = v.planArea() * (v.hi.z - v.interfaceZ);
        expectNear("and the drawn hot-layer volume is the model's upper volume", drawnUpper,
                   g.upperVolume(), std::max(1e-9, g.upperVolume() * 1e-9));
        expectTrue("the interface sits inside the space",
                   v.interfaceZ >= v.lo.z - 1e-12 && v.interfaceZ <= v.hi.z + 1e-12);

        // **The plan area alone does not pin the prism.** A square footprint of the
        // same area satisfies every identity above -- the volume, the layer split
        // and the interface are all right -- and draws a 24 m engine room as a 15 m
        // square. The aspect ratio is the compartment's, and it is checked here
        // because nothing else looks at it.
        const sim::Compartment& c =
            ship.compartments[static_cast<std::size_t>(g.shipCompartment)];
        const double wanted = (c.bboxHi.x - c.bboxLo.x) / (c.bboxHi.y - c.bboxLo.y);
        expectNear("and the prism keeps the compartment's plan aspect ratio",
                   (v.hi.x - v.lo.x) / (v.hi.y - v.lo.y), wanted, wanted * 1e-12);
        expectNear("centred on the compartment, in x",
                   0.5 * (v.lo.x + v.hi.x), 0.5 * (c.bboxLo.x + c.bboxHi.x), 1e-9);
        expectNear("and in y", 0.5 * (v.lo.y + v.hi.y), 0.5 * (c.bboxLo.y + c.bboxHi.y), 1e-9);
    }
    // Vacuity: the ferry's engine rooms are nothing like square, so the aspect
    // check above is a real constraint rather than one a square would satisfy.
    const gpu::SmokeVolume& room = volumes[static_cast<std::size_t>(
        model.findGas("engine_room_s"))];
    expectTrue("and the compartment it was checked on is not square",
               (room.hi.x - room.lo.x) > 2.0 * (room.hi.y - room.lo.y));
}

// The lower layer's optics come from the lower layer's own volume, on the path a
// caller actually takes.
//
// **The ferry cannot ask this either.** Its cool layer never carries any tracer,
// so `volumesFromFire` reading the whole gas volume for it -- rather than what is
// left after the hot layer -- produces exactly the same zero. The state is set by
// hand here for the same reason `Model::gas` is public: a model assembled for a
// question is the way to ask one the ship does not.
void testTheCoolLayerIsSizedByWhatIsLeft() {
    sim::Ship ship = game::buildFerry();
    const sim::Sea sea;
    ship.initialise(sea);
    sim::fire::Model model = ferryFire(ship, sea, 300);
    const int burning = model.findGas("engine_room_s");
    sim::fire::GasCompartment& g = model.gas[static_cast<std::size_t>(burning)];

    // Enough tracer in the cool layer to be seen, and a different amount from the
    // hot one so the two cannot be confused.
    g.lower.products = 0.25;
    const double upperVolume = g.upperVolume();
    const double lowerVolume = g.gasVolume - upperVolume;
    expectTrue("the fixture has two layers of different size",
               upperVolume > 1.0 && lowerVolume > 1.0 &&
                   std::abs(upperVolume - lowerVolume) > 0.2 * g.gasVolume);

    const gpu::SmokeShading shading;
    const std::vector<gpu::SmokeVolume> volumes = gpu::volumesFromFire(model, ship, shading);
    const gpu::SmokeLayer& lower = volumes[static_cast<std::size_t>(burning)].lower;
    expectNear("the cool layer's extinction is its own products over its own volume",
               lower.extinction, shading.massExtinction * g.lower.products / lowerVolume,
               lower.extinction * 1e-12);
    // And it is *not* what the whole gas volume would give, by a wide margin.
    expectTrue("which is not what the whole gas space would give",
               std::abs(lower.extinction -
                        shading.massExtinction * g.lower.products / g.gasVolume) >
                   0.1 * lower.extinction);
}

// The optical chain, end to end, on a real fire: concentration from the model's
// own layer volume, extinction from the tracer, emission from the temperature.
void testTheOpticsFollowTheGasState() {
    sim::Ship ship = game::buildFerry();
    const sim::Sea sea;
    ship.initialise(sea);
    const sim::fire::Model model = ferryFire(ship, sea, 400);
    gpu::SmokeShading shading;
    const std::vector<gpu::SmokeVolume> volumes = gpu::volumesFromFire(model, ship, shading);

    const int burning = model.findGas("engine_room_s");
    const int next = model.findGas("engine_room_p");
    expectTrue("the ferry fixture has both engine rooms", burning >= 0 && next >= 0);
    const gpu::SmokeVolume& fire = volumes[static_cast<std::size_t>(burning)];
    const gpu::SmokeVolume& door = volumes[static_cast<std::size_t>(next)];
    const sim::fire::GasCompartment& g = model.gas[static_cast<std::size_t>(burning)];

    // Extinction is the mass extinction coefficient times a concentration, and
    // the concentration is the model's, not the box's.
    const double concentration = g.upper.products / g.upperVolume();
    expectNear("the hot layer's extinction is K_m times its own concentration",
               fire.upper.extinction, shading.massExtinction * concentration,
               fire.upper.extinction * 1e-12);
    expectTrue("which is an opaque layer on this fire", fire.upper.extinction > 1.0);
    expectTrue("the cool layer is far clearer",
               fire.lower.extinction < 0.2 * fire.upper.extinction);

    // Species crossed the door, so the room next to the fire is smoky and dimmer.
    expectTrue("smoke reached the compartment next door", door.upper.extinction > 1e-3);
    expectTrue("but it is thinner there", door.upper.extinction < fire.upper.extinction);

    // Emission follows temperature, and the burning room is the hot one.
    double predicted[3];
    gpu::emissiveColour(g.upper.temperature(), shading, predicted);
    for (int c = 0; c < 3; ++c)
        expectNear("the hot layer emits the Planck colour of its own temperature",
                   fire.upper.emission[c], predicted[c],
                   std::max(1e-18, predicted[c] * 1e-12));
    expectTrue("the burning room's layer is hotter than its neighbour's",
               fire.upper.emission[0] > door.upper.emission[0]);
    // Vacuity in the other direction: the fire has to have made the layer hot
    // enough for the check above to be about anything.
    expectTrue("and the fire really did heat it", g.upper.temperature() > 500.0);

    // The knob is a knob: doubling K_m doubles every extinction and touches no
    // emission. A coefficient that quietly scaled the colour too would be a
    // second, hidden exposure control.
    shading.massExtinction *= 2.0;
    const std::vector<gpu::SmokeVolume> doubled = gpu::volumesFromFire(model, ship, shading);
    expectNear("doubling the mass extinction coefficient doubles the extinction",
               doubled[static_cast<std::size_t>(burning)].upper.extinction,
               2.0 * fire.upper.extinction, fire.upper.extinction * 1e-12);
    expectNear("and leaves the emission alone",
               doubled[static_cast<std::size_t>(burning)].upper.emission[0],
               fire.upper.emission[0], 1e-15);
}

// The layer descends, and the picture has to descend with it. A renderer that
// read the interface once and cached it would pass every static check above.
void testTheDrawnInterfaceFollowsTheLayerDown() {
    sim::Ship ship = game::buildFerry();
    const sim::Sea sea;
    ship.initialise(sea);
    sim::fire::Model model;
    model.attach(ship, {ship.findCompartment("engine_room_s"),
                        ship.findCompartment("engine_room_p")});
    sim::fire::DesignFire d;
    d.name = "machinery";
    d.compartment = model.findGas("engine_room_s");
    d.baseZ = 2.5;
    d.diameter = 2.5;
    d.growthCoefficient = sim::fire::kGrowthFast;
    d.peakHeatRelease = 4.0e6;
    d.steadyDuration = 900.0;
    model.fires.push_back(d);

    const gpu::SmokeShading shading;
    const int burning = model.findGas("engine_room_s");
    std::vector<double> heights, opacity;
    for (int block = 0; block < 12; ++block) {
        for (int i = 0; i < 50; ++i) model.step(1.0, ship, sea);
        const std::vector<gpu::SmokeVolume> v = gpu::volumesFromFire(model, ship, shading);
        heights.push_back(v[static_cast<std::size_t>(burning)].interfaceZ);
        opacity.push_back(v[static_cast<std::size_t>(burning)].upper.extinction);
    }

    // **The layer does not descend for ever, and the first version of this test
    // asserted that it did.** On the ferry it reaches 2.95 m at about t = 300 s
    // and then *recovers* to 3.13 m as the room reaches its vented steady state
    // and the hot layer's mass falls. That is the model's answer and the renderer
    // has to follow it, so what is asserted is monotone descent through the growth
    // phase and a bounded recovery after it -- not a monotonicity the physics does
    // not have.
    // Index 5 is t = 300 s, the last sample before the recovery: measured, by
    // sampling this run every ten seconds and reading where the minimum is.
    const std::size_t growth = 5;
    bool descending = true, thickening = true;
    for (std::size_t i = 1; i < heights.size(); ++i) {
        if (i <= growth && heights[i] > heights[i - 1] + 1e-9) descending = false;
        if (opacity[i] < opacity[i - 1]) thickening = false;
    }
    expectTrue("the drawn interface descends throughout the growth phase", descending);
    expectTrue("and the hot layer only ever gets more opaque, over the whole run", thickening);
    expectTrue("and it moved far enough to see", heights.front() - heights[growth] > 0.5);
    expectTrue("the recovery afterwards is small", heights.back() < heights[growth] + 0.5);
    // Vacuity: an extinction that never left zero would satisfy "only ever
    // increases" perfectly.
    expectTrue("and the layer became genuinely opaque", opacity.back() > 1.0);
}

}  // namespace

void runSmokeTests() {
    std::printf("\n--- volumetric fire and smoke: the optics ---\n");
    testTheSpectrumIntegratesToStefanBoltzmann();
    testPlanckReachesItsRayleighJeansLimit();
    testThePeakSitsWhereWienSaysItDoes();
    testTheBandsAreASpectralSampleAndTheMappingIsOneNumber();
    testTransmittanceIsBeerLambert();
    testTheRayMeetsThePrismWhereGeometrySaysItDoes();
    testThePathSplitsAtTheInterface();
    testTheMediumRotatesWithTheShip();
    testTheCompositeIsTheTransferEquation();
    testTheDepthBasisComesOutOfTheMatrix();
    testBothLayersTakeTheirOwnConcentration();
    testTheDrawnVolumeIsTheModelsVolume();
    testTheCoolLayerIsSizedByWhatIsLeft();
    testTheOpticsFollowTheGasState();
    testTheDrawnInterfaceFollowsTheLayerDown();
}
