// SPDX-License-Identifier: MIT
//
// Tier 1: Craig-Bampton component mode synthesis -- the missing middle of the
// structural model.
//
// `girder.hpp` is the whole ship as a beam: 200 DOF, microseconds, and it knows
// nothing about where stress goes *within* a section. `zone.hpp` is solid-shell
// elements with plasticity over a patch: uncompromised, and `4.0 x elementCount`
// core-seconds per simulated second with a hard ceiling at ~8000 elements, where
// the per-element stiffness store leaves L3 and the parallel speedup falls off a
// cliff (`02-simulation.md` §3). Between "the ship is a beam" and "twenty-four
// square metres of plating" there was nothing, and a whole hold, a whole
// superstructure, or the region between two bulkheads is exactly the size that
// falls in the gap.
//
// This file closes it. A substructure's degrees of freedom are partitioned into
// **boundary** (the interface it shares with the rest of the ship) and
// **interior**. Every boundary DOF is kept *exactly*. The interior -- which is
// almost all of them -- is represented by the static constraint modes plus a
// handful of fixed-interface normal modes. What comes out is a small dense mass
// and stiffness pair that reproduces the substructure's response at its interface
// to whatever accuracy the mode count buys, and that a caller can integrate at a
// step set by the physics rather than by the explicit stability limit.
//
// --- 1. The formulation --------------------------------------------------------
//
// Partition u = [u_i ; u_b] and write
//
//     u_i = Psi u_b + Phi q,      Psi = -K_ii^-1 K_ib
//
// `Psi` is the **constraint modes**: the exact static interior displacement for a
// unit displacement of each interface DOF in turn. `Phi` is the **fixed-interface
// normal modes**: the lowest `m` eigenvectors of (K_ii, M_ii), mass-normalised, so
// Phi^T M_ii Phi = I and Phi^T K_ii Phi = Lambda = diag(omega_j^2). The
// transformation is
//
//     u = T x,   x = [u_b ; q],   T = [[Psi, Phi], [I, 0]]
//
// and the reduced pair is K_r = T^T K T, M_r = T^T M T:
//
//     K_r = [[K_bb - K_bi K_ii^-1 K_ib,  0     ],   <- the Schur complement; the
//            [0,                         Lambda]]      off-diagonal is *identically*
//                                                      zero, which is what makes
//                                                      the reduction stable.
//     M_r = [[M_bb + Psi^T M_ii Psi,     Psi^T M_ii Phi],
//            [Phi^T M_ii Psi,            Phi^T M_ii Phi]]
//
// **Two of the stiffness blocks are stored as their closed forms rather than
// computed**, and one of the mass blocks is not. The (b,q) stiffness block being
// zero is not an approximation: K_bi Phi + Psi^T K_ii Phi = K_bi Phi
// - K_bi K_ii^-1 K_ii Phi = 0, so it is stored as the zero it is; and the modal
// block is Lambda, the eigenvalue, where computing Phi^T K_ii Phi instead would
// be a noisier route to the same number. The modal *mass* block is the identity
// by the mass normalisation, and that one **is** computed, because it is cheap and
// because asserting on it is what says the eigenvectors really came back
// mass-normalised. `tests/test_reduction.cpp` checks all of it by forming T^T K T
// and T^T M T the long way and comparing.
//
// **Four properties are exact and are tested as identities, not tolerances:**
//
//  1. **Zero modes is Guyan static condensation, and static condensation is exact
//     at the interface.** For *any* load -- not just one applied at the boundary --
//     u_b = Khat^-1 (f_b + Psi^T f_i) is the same solve as eliminating the interior
//     from the full system. The mode count buys nothing at the interface
//     statically; what it buys is the *interior* recovery under an interior load,
//     and dynamics. (§5 has the measurement, and this corrects the usual statement
//     of the property, which is that the static interface response "improves with
//     mode count". It does not: it starts exact.)
//  2. **Rigid body modes survive exactly.** A rigid motion of the interface is the
//     exact static solution of the interior with no interior load, so
//     `Psi u_b^rigid` *is* the rigid interior field and Khat u_b^rigid = 0. A
//     free-free substructure keeps exactly six zero eigenvalues -- measured at
//     2e-10 of the first elastic one -- and because the mass goes through the same
//     T, x^T M_r x for a rigid translation is the substructure's own total mass and
//     for a rigid rotation is the full model's own rotary inertia.
//  3. **Reduced frequencies come down from above.** The reduction is a
//     Rayleigh-Ritz projection onto a subspace of the full space, so every reduced
//     eigenvalue is an upper bound on the corresponding full one, and the
//     subspaces nest as modes are added, so the bound falls monotonically. A
//     reduction can only stiffen. The tests assert the *direction*, because
//     "close" would pass for a reduction that had the sign of its error wrong.
//  4. **Symmetry and definiteness carry over -- but not the same definiteness for
//     both matrices, and the difference is property 2.** M_r is positive
//     *definite*: T has full column rank and M is a positive diagonal, so
//     T^T M T is too. K_r is positive **semi**-definite, and for a free-free
//     substructure it is exactly six-fold singular -- being positive definite
//     would mean the rigid body modes had been lost. Both are symmetric to the
//     last bit, because the analytically symmetric blocks are computed on one
//     triangle and mirrored; that makes the symmetry *check* vacuous on its own,
//     which is why the test that carries the content is the independent T^T K T.
//
// --- 2. Why an eigensolver had to be written -----------------------------------
//
// No third-party dependencies, the same reason this repo has a hand-written PNG
// codec and its own Threefry. Two solvers, because the two problems are different
// shapes:
//
//   * **Dense**, for the projected q x q problem inside the subspace iteration and
//     for anything a test wants the whole spectrum of: Householder reduction to
//     tridiagonal form with the transform accumulated, then implicit QL with
//     Wilkinson shifts. O(n^3), unconditionally convergent, and the standard
//     answer.
//   * **Subspace iteration**, for the fixed-interface modes themselves, because
//     only the lowest `m` of several thousand are wanted and a dense solve of the
//     interior block is O(n_i^3). It iterates X <- K_ii^-1 M_ii X on a block of
//     `q > m` vectors and solves the small projected problem with the dense solver
//     above. The banded Cholesky of K_ii is `solidshell::BandedSpd`, which already
//     exists and is already validated.
//
// **The block is re-orthogonalised in the mass inner product every iteration, and
// that is not a refinement.** Without it, repeated multiplication by K^-1 M drives
// every column towards the lowest modes and the block goes numerically rank
// deficient: measured, on a 210-DOF interior with a block of 200, the projected
// mass matrix stopped being positive definite, the projected eigenproblem was
// refused, and *the whole reduction fell back to Guyan* with one line in
// `problems` to say so. With the block orthonormal the projected mass matrix is the
// identity by construction and there is no Cholesky left to fail -- so the
// projected problem is a standard symmetric one rather than a generalised one. It
// is not free: it is what takes a 47-mode ferry patch from 0.33 s to 1.5 s.
//
// **Subspace iteration converges to the lowest modes; it does not prove it found
// them.** A skipped mode is the silent failure here -- the reduced model would
// still be symmetric, still positive definite, still bound the full model from
// above, and simply be missing a mode. So the mode count is verified by a
// **Sturm sequence check**: an LDL^T factorisation of (K_ii - sigma M_ii) has
// exactly as many negative pivots as there are eigenvalues below sigma
// (Sylvester's law of inertia). Put sigma between the last mode kept and the first
// discarded and the count must equal `m` exactly. `Reduction::modesVerified` says
// whether it did.
//
// The same inertia count is what makes a **frequency cutoff** an exact request
// rather than a guess: "every mode below 20 Hz" is answered by one factorisation
// before any eigenvector is computed.
//
// --- 3. How the interface is chosen ---------------------------------------------
//
// **The interface is every node the substructure shares with anything outside it,
// and nothing else.** That is a property of the cut, not of the mesh, so the
// caller supplies it; a mesh handed over in isolation has no "outside" to infer
// one from. Two helpers cover the cases that actually arise -- `nodesNearPlanes`
// for a region cut out between two bulkheads, and `nodesPinned` for a
// `zone::buildPatch` patch, whose clamped perimeter is exactly the set of nodes
// the surrounding plating holds.
//
// Two rules the caller has to respect, and both are checked:
//
//   * **The interface must restrain the interior**, and the factorisation does
//     *not* catch it when it does not. K_ii is the stiffness with the interface
//     held, so it is singular unless the interface pins all six rigid body modes.
//     A mechanism leaves an exactly zero pivot in exact arithmetic and a tiny
//     **positive** one in floating point, so `BandedSpd::factor` succeeds and the
//     solve returns nonsense -- measured, on a two-node interface. The
//     precondition is geometric, so it is checked geometrically before anything is
//     factored: at least three nodes not on a line, with the tolerance scaled by
//     the interface's own extent so that one collinear to within a rounding is
//     refused too. (A mechanism *inside* the mesh would survive it, and a caller
//     who suspects one can count `eigenvaluesBelow` a small positive shift.)
//   * **Cost is quadratic in the interface, linear in the interior.** Psi is
//     n_i x n_b dense and Psi^T M_ii Psi is O(n_i n_b^2), so an interface twice as
//     big costs four times as much and stores twice as much. A bulkhead cut is
//     cheap because a transverse section is a thin ring of nodes; "every node on
//     the outside of the hold" would not be, and would also be the wrong cut.
//
// **`HexMesh::fixed` is ignored, deliberately.** A substructure is a free
// component; what holds it is a constraint on the *reduced* model, applied by the
// caller through `staticSolve`'s held-DOF list. Silently consuming the mesh's own
// pins would make the free-free rigid body property untestable and would mean a
// patch could not be reduced with its clamped edge as its interface -- which is
// the main thing anyone wants to do with one. `Substructure::problems()` says so
// when the mesh carries pins, because a dropped constraint that is not reported is
// indistinguishable from one that was honoured.
//
// --- 4. Mass is reduced too, through the same T ---------------------------------
//
// Yes, and there is no defensible alternative. The two things that could be done
// instead are both measurably wrong:
//
//   * **Keeping M_bb alone** (the "Guyan mass" a static condensation leaves) hands
//     the reduced model only the interface's share of the inertia. On the test
//     substructures here that is 20-30% of the ship it stands for, so a rigid
//     translation of the reduced model reports a fraction of the real mass, and
//     every frequency comes out high by the square root of the ratio.
//   * **Reducing K through T and M through anything else** loses the bound in
//     property 3 above: the from-above convergence is a Rayleigh quotient
//     statement about *one* subspace, and it is only true when both matrices are
//     projected onto it.
//
// (§8 adds one term to this: a stiffener condensed onto the plating's degrees of
// freedom brings mass that no element carries, and it is lumped onto the same
// diagonal for the same reasons.)
//
// The mass matrix itself is **row-sum lumped**, matching Tier 2, and that is a
// choice with consequences worth stating. In its favour: `zone.hpp`'s explicit
// solver lumps, so the reduced region and the zone that replaces it agree about
// inertia rather than disagreeing at the moment of promotion; a diagonal M_ii
// turns the fixed-interface eigenproblem into a standard symmetric one with no
// Cholesky of a mass matrix; M_bi vanishes identically, which removes a whole
// block from the reduction; and row-sum lumping conserves total mass and rigid
// translation exactly. Against it: a lumped mass under-states rotary inertia, so
// the *absolute* frequencies sit a little below a consistent-mass model's. That
// does not touch anything asserted here, because every comparison is against the
// **same** full finite element model, and the from-above property is a statement
// about the discrete pair, not about the continuum.
//
// --- 5. How many modes, and what they buy --------------------------------------
//
// The rule is a **frequency cutoff at twice the highest frequency of interest**,
// which is the standard component-mode-synthesis practice and is what
// `ReduceParams::cutoffFrequency` defaults to. The band this engine actually has
// loads for is set by the hull girder: a 120 m ferry's two-node vertical mode is
// around 1.5 Hz and the fifth is around 10 Hz, whipping and springing live in
// that band, so 10 Hz is the frequency of interest and 20 Hz (126 rad/s) is the
// cutoff. A slam transient has content an order higher, and a caller who wants one
// resolved must say so -- there is no default that is right for both.
//
// **What the rule actually buys, measured, is about one per cent -- not the "few
// per cent at 10^-5 of the cost" `02-simulation.md` §3 first claimed and not the
// four figures the smoothness of the convergence suggests.** On a 272-element
// patch of the ferry's own side (1902 DOF, 528 of them interface) with half the
// interface free, against a dense solve of the same full model:
//
//     modes  fixed-interface   worst error over    worst error over
//            cutoff it is      the modes < 10 Hz   the modes < 20 Hz
//     ------------------------------------------------------------------
//        0   (Guyan alone)          3.9e-2              4.9e-1
//        2        16 Hz             6.0e-3              8.3e-2
//        5        24 Hz             5.0e-3              1.5e-2
//       10        34 Hz             3.8e-3              5.3e-3
//       20        56 Hz             5.7e-4              1.9e-3
//       47       116 Hz             8.8e-5              2.0e-4
//      100       301 Hz             3.7e-6              8.3e-6
//
// The trap in reading that is the one this file walked into: **the cutoff is a
// frequency of the fixed-interface problem, and the band of interest is a
// frequency of the assembled one, and they are not the same spectrum.** With part
// of the interface free the assembled model is much softer than the
// fixed-interface one -- 2.97 Hz against a first fixed-interface mode at 14.7 Hz --
// so "everything below 2 x 10 Hz" keeps *two* modes here, not the dozens the
// frequency count of the assembled model would suggest. Two modes is 0.6% inside
// the 10 Hz band and 8% up to 20 Hz. Five or six times the band, rather than
// twice, is what buys a part in ten thousand.
//
// The second measured consequence determines when this file is worth calling at
// all: **a small stiff substructure needs no modes.** Its first fixed-interface
// frequency is already far above the band, so Guyan condensation is the whole
// answer -- and it is exactly right statically at the interface however soft the
// load. `Reduction::firstFixedFrequency` is computed and reported even when zero
// modes are kept, precisely so that "no modes were needed" is a number rather than
// an assumption.
//
// --- 6. What a reduced region cannot do ----------------------------------------
//
// **A reduced model is linear by construction, and that is not a limitation that
// can be worked around inside it.** K is formed once from the undeformed
// geometry; Psi and Phi are a fixed subspace derived from that K. Consequently:
//
//   * **It cannot yield.** There is no stress state to return-map, no plastic
//     history, no path dependence. Feed it a load past yield and it will report a
//     linear-elastic answer that is too stiff and too strong, without complaint.
//   * **It cannot tear, buckle, or contact anything.** Element deletion, geometric
//     stiffening, and the whole of `plasticity.hpp` are outside the subspace.
//   * **It cannot rotate.** The stiffness is small-strain in the global frame, so
//     a substructure that turns through a finite angle -- a capsizing ship, a
//     detached section -- is outside its validity, unlike the co-rotational Tier-2
//     element.
//   * **`checkValidity` reads the elements only, so an attached member's own
//     stress is not in the utilisation** -- see §8. A stiffened region is reported
//     by its plating, and the plating is the softer half.
//   * **It is stale the moment anything changes.** Thinning, tearing or a moved
//     boundary changes K, so the reduction must be rebuilt. It is a one-off cost,
//     not a per-step one, but it is not free.
//
// **What the caller must do instead.** Watch `checkValidity`, which recovers the
// interior displacement, evaluates the element stresses through the same
// `solidshell::elementStress` the Tier-2 solver uses, and reports the peak von
// Mises utilisation. When it approaches one, the region has to be **promoted to
// Tier 2** -- `promotion.hpp` already owns that decision and `zone::buildPatch`
// already meshes the result. The reduced model's job at that point is to stop
// answering and hand over its interface displacements as the zone's boundary
// condition; it must not be asked for the answer itself. **Handing them over is
// `coupling.{hpp,cpp}`**, which drives a `zone::Patch`'s perimeter from an
// assembled reduced model and is exact for a linear zone -- the property below
// being exactly why.
//
// One caveat on that warning, and it points the wrong way: the recovered stress is
// only as good as the mode set. A truncated basis cannot represent a stress
// concentration, so the peak comes out **under**-predicted and the warning is
// late rather than early. `checkValidity` is a trigger for promotion, not a
// strength check, and the measurement of how late it is is in
// `02-simulation.md` §3.
//
// --- 7. What it costs, against the tiers either side of it ----------------------
//
// Measured on the same 272-element patch of the ferry's side, so the three numbers
// are comparable rather than merely quoted:
//
//   * **Tier 0**, the whole ship as a beam: `hullGirder` is 0.995 ms, so 0.10
//     core-seconds per simulated second at the 100 Hz §3 specifies -- for the
//     *whole ship*, not for this patch.
//   * **Tier 1**, this patch: 0.11 s to reduce with no modes, 1.5 s with 47, plus
//     0.01 s to factor the implicit operator once. Then **0.35-0.41 core-seconds
//     per simulated second** at 1 ms steps, because a linear model is
//     unconditionally stable and its step is set by what has to be resolved.
//   * **Tier 2**, the same patch: **1155 core-seconds per simulated second**. The
//     explicit step is 1.72 us -- thickness-governed, not chosen -- so it takes
//     580 000 steps to Tier 1's thousand.
//
// So Tier 1 is **2800x cheaper than Tier 2 on the same plating** and about four
// times the cost of the beam that covers the entire ship. That is the trade: three
// orders of magnitude of speed for everything nonlinear.
//
// **The interface sets the running cost, not the mode count.** Going from 0 to 47
// modes moves a step from 350 to 409 us -- 17% -- while the 528 interface degrees
// of freedom are 92% of the reduced model either way. An interface twice as large
// costs four times as much to build and twice as much to store, which is the
// argument for cutting a substructure at a bulkhead rather than around a patch.
//
// --- 8. Stiffness and mass the elements do not carry ----------------------------
//
// A `HexMesh` is elements, and **not everything structural is an element**. The
// case this exists for is the one `constraint.hpp` builds: an eccentric stiffener
// is represented as axial fibres tied to the plating and *condensed onto the
// plating's own degrees of freedom*, so it has no nodes and no elements of its
// own. A substructure that only ever looked at the mesh therefore reduced a
// stiffened patch as **bare plating** -- measured, before this section existed, at
// 3e-16 relative to an element-only assembly, while one 200 x 10 flat bar across a
// 0.6 m patch was worth 7.9% of its displacement field. Longitudinals carry a
// large share of a hull girder, so that error is not a patch-scale curiosity: it
// is wrong in the unsafe direction for any Tier-1 model of real ship structure.
//
// `Attachment` closes it. It is deliberately expressed in `solidshell::DofBlock`
// and a nodal mass array rather than in `constraint::Stiffening`, for two reasons:
// the substructure has no business knowing what a stiffener is, and `DofBlock` is
// already exactly what `solidshell::solveStatic` takes -- so the static reference
// a test compares against is the *same* description of the *same* physics, not a
// second one that could disagree. `constraint::stiffnessBlocks` and
// `constraint::lumpFiberMass` produce the two fields between them.
//
// **Three things had to be got right, and each is a silent failure otherwise.**
//
//   1. **The sparsity pattern.** The CSR is built from element node adjacency, and
//      a block ties degrees of freedom that no single element need share. The
//      lookup is a binary search that, on a column the pattern does not have,
//      returns the slot of the *next* one -- so an unchecked scatter would land a
//      stiffener's stiffness on a neighbouring degree of freedom and produce a
//      plausible field that is quietly wrong. The blocks are therefore folded into
//      the adjacency *before* the pattern is built, and every scatter checks that
//      the slot it found is the slot it asked for.
//   2. **The bandwidth.** `solidshell::BandedSpd::add` drops anything outside the
//      band it was given, without a word. The band here is taken from the
//      assembled pattern rather than estimated from the elements, so folding the
//      blocks in above fixes this too -- and the interior renumbering sees them,
//      which it would not if the blocks were scattered after the fact. That it
//      worked is checked rather than argued: the banded fill counts what fell
//      outside the band and refuses the substructure if anything did.
//
//      **How far outside the band a block actually reaches is worth knowing,
//      because the obvious answer is wrong.** A fibre from
//      `constraint::addStiffener` spans two *consecutive* stations of a seam, and
//      consecutive stations are corners of the same element -- so at mesh
//      resolution every entry a stiffener contributes is already in the pattern
//      and the band does not move at all (measured: 41 to 41 on the test plate).
//      It moves when the member is modelled *coarser* than the plating it sits on:
//      a station every second node takes the same plate 41 to 65. That is the case
//      the tests use, because a test built only on the default meshing would have
//      proved nothing about any of this.
//   3. **Mass.** A block is stiffness alone; the steel it stands for weighs
//      something. Leaving it out gives a model that is stiffer *and no heavier*,
//      so every frequency comes out high -- and frequencies are most of what this
//      tier is for. So the mass is taken as well, lumped onto the same diagonal as
//      the elements' (§4), which is also what `zone::Solver` does with
//      `constraint::lumpFiberMass` -- the two therefore agree about inertia at the
//      moment of promotion rather than disagreeing. A substructure given stiffness
//      and no mass says so in `problems()`.
//
// **What is given up, and it is `constraint.hpp`'s measurement rather than a new
// one.** The eccentric tie's consistent condensed mass `T^T M T` puts a *negative*
// mass on one of the two through-thickness nodes for any weight outside [0, 1] --
// -7.83 m against +8.83 m for a 200 mm bar on 12 mm plating -- so the fibre mass is
// split equally over the pair instead. The total is exact and the rotary inertia
// about the seam is not represented. Tier 1 is implicit and could integrate a
// negative diagonal entry where Tier 2 could not, but taking the consistent mass
// here and the lumped one there would make the two tiers disagree about the
// inertia of the same steel, which is the disagreement §4 exists to avoid.
//
// **And what it does not reach.** `checkValidity` walks the elements, so an
// attached member's own stress is not in the peak von Mises -- a stiffened region
// is judged by its plating. The stiffener makes the plating *less* utilised, so
// the promotion trigger moves later on both counts, which is the same direction
// §6 already warns about for a truncated basis. Measured on the stiffened
// cantilever in the tests: 58.2 MPa reported against 65.2 MPa in the member, a
// utilisation of 0.164 where the true one is 0.184.
//
// **One result from validating this is worth carrying, because it is not the
// expected one.** A stiffener stiff enough holds its own line still, and that line
// then becomes a *node* of the first mode: on the test plate the 200 x 10 bar
// raises the first frequency by 307% while its own mass moves that same frequency
// by 0.015%. The mass is not negligible -- 13.7% on the worst of the first twelve
// modes -- but a single-frequency check would have concluded it did not matter,
// and on a member the panel can actually bend with (60 x 6) the first frequency
// lands within 4.5% of what `scantlings::stiffenedSection` predicts.
//
// SI units, body frame per CLAUDE.md.
#pragma once

#include "scantlings.hpp"
#include "solid_shell.hpp"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace sim::reduction {

// --- Symmetric eigensolvers ----------------------------------------------------

// Eigenpairs, ascending by eigenvalue. `vector` holds `size` doubles per pair,
// pair j at [j * size, (j + 1) * size).
struct Eigenpairs {
    std::vector<double> value;
    std::vector<double> vector;
    int size = 0;        // length of each eigenvector
    int count = 0;       // how many pairs
    int iterations = 0;  // QL sweeps, or subspace iterations
    bool converged = false;
    std::string problem;

    const double* mode(int j) const { return vector.data() + static_cast<std::size_t>(j) * static_cast<std::size_t>(size); }
};

// A x = lambda x for symmetric A, `n` x `n` row-major. Only the lower triangle is
// read, so a caller with a symmetric matrix stored either way gets the same
// answer. Householder tridiagonalisation then implicit QL with Wilkinson shifts;
// eigenvectors are orthonormal.
Eigenpairs symmetricEigen(const std::vector<double>& a, int n);

// A x = lambda B x with B symmetric positive definite. Reduced to standard form
// through the Cholesky factor of B, so the eigenvectors come back **B-orthonormal**
// -- x_i^T B x_j = delta_ij -- which is the normalisation a modal model needs.
// `converged` is false with a reason in `problem` when B is not positive definite.
Eigenpairs generalisedEigen(const std::vector<double>& a, const std::vector<double>& b, int n);

// --- Choosing an interface ------------------------------------------------------

struct Plane {
    Vec3 point{};
    Vec3 normal{};  // need not be unit; it is normalised here
};

// Nodes within `tolerance` of any of the planes: the region between two bulkheads,
// cut at the two transverse sections it is attached through. Ascending, unique.
std::vector<std::uint32_t> nodesNearPlanes(const solidshell::HexMesh& mesh,
                                           const std::vector<Plane>& planes,
                                           double tolerance = 1e-9);

// A bandwidth-reducing permutation of a node adjacency graph: reverse Cuthill-McKee
// from a George-Liu pseudo-peripheral start. `order[i]` is the node that belongs at
// position `i`, and every node appears once.
//
// It is used here on the interior partition, and it is exposed because
// `section.hpp` needs the same thing for a different reason --
// `solidshell::solveStatic` numbers its free degrees of freedom in the **mesh's own
// order** and has no renumbering pass at all, so a section mesher has to arrive with
// a good one or pay a hundredfold in the factorisation. A second copy of this would
// be a second place the two could disagree.
//
// **It is not unconditionally better than the numbering it replaces**, which is
// measured rather than assumed: 59 against 41 on one test mesh and 83 against 53 on
// a ferry patch, against 89 to 173 and 137 to 341 on two others. Both callers
// therefore compare it against what they already had and keep the narrower, which is
// free -- a bandwidth is one pass over the adjacency.
std::vector<std::uint32_t> bandwidthReducingOrder(
    const std::vector<std::vector<std::uint32_t>>& adjacency);

// Nodes any of whose DOF the mesh pins. A `zone::buildPatch` patch arrives with
// its perimeter clamped, and that perimeter is exactly the set of nodes the
// plating outside the patch holds -- so this turns a Tier-2 patch into a
// substructure whose interface is the boundary the patch already has.
std::vector<std::uint32_t> nodesPinned(const solidshell::HexMesh& mesh);

// --- The substructure ----------------------------------------------------------

// Stiffness and mass that belong to the component but not to any of its elements:
// a stiffener condensed onto the plating it is welded to, and the steel that
// stiffener is made of. See §8 for what this is for and what it costs.
//
// Both fields come ready-made from `constraint.{hpp,cpp}`:
//
//     const constraint::RestFibers forms = constraint::restFibers(stiffening, mesh.position);
//     reduction::Attachment attached;
//     attached.stiffness = constraint::stiffnessBlocks(stiffening, mesh.position, forms,
//                                                      material.youngsModulus);
//     attached.mass.assign(mesh.nodeCount(), 0.0);
//     constraint::lumpFiberMass(stiffening, forms, material.density, attached.mass);
//
// **`mass` is per node, not per degree of freedom**, which is the convention
// `constraint::lumpFiberMass` and `zone::Solver` already use; a wrongly sized
// array is refused rather than read short, because an array three times too long
// would otherwise assemble a plausible model with two thirds of the steel missing.
// Empty means no attached mass, which is a defensible choice for a pure stiffness
// tie and is reported so it cannot be an oversight.
struct Attachment {
    std::vector<solidshell::DofBlock> stiffness;
    std::vector<double> mass;  // kg per node, added to the row-sum lumped diagonal

    bool empty() const { return stiffness.empty() && mass.empty(); }
};

// The full finite element model of one component, assembled once, partitioned, and
// factored on its interior. Holds a **reference** to the mesh: it must outlive the
// substructure, in the same way `zone::Solver` holds its `Patch`.
class Substructure {
public:
    // Elements alone. Delegates to the constructor below with an empty
    // `Attachment`, so there is exactly one assembly path and "no attachment
    // reduces to what this always did" is a property of the code rather than a
    // promise a test has to keep checking -- though `tests/test_reduction.cpp`
    // checks it anyway, bit for bit, because that is the negative control for
    // everything §8 adds.
    Substructure(const solidshell::HexMesh& mesh, const StructuralMaterial& material,
                 const std::vector<std::uint32_t>& interfaceNodes,
                 solidshell::Formulation form = solidshell::Formulation::SolidShell);
    // Elements plus whatever else the component carries. See §8. The `Attachment`
    // is **consumed here, not held** -- unlike the mesh, which is a reference --
    // because its blocks are scattered into the same store the elements go into and
    // there is nothing left to refer back to. It therefore need not outlive the
    // call, and changing it afterwards changes nothing.
    Substructure(const solidshell::HexMesh& mesh, const StructuralMaterial& material,
                 const std::vector<std::uint32_t>& interfaceNodes, const Attachment& attached,
                 solidshell::Formulation form = solidshell::Formulation::SolidShell);
    ~Substructure();
    Substructure(Substructure&&) noexcept;
    Substructure& operator=(Substructure&&) noexcept;

    // False when the interior could not be factored -- almost always an interface
    // that does not restrain all six rigid body modes -- or the mesh is degenerate.
    bool ready() const;
    const std::vector<std::string>& problems() const;

    std::size_t nodeCount() const;
    std::size_t dofCount() const;      // 3 per node
    std::size_t boundaryCount() const;
    std::size_t interiorCount() const;
    std::size_t halfBandwidth() const;  // of K_ii after the interior renumbering
    double assemblySeconds() const;

    // What the `Attachment` contributed. `attachedMass` is the sum of the array
    // that was handed over, so a caller can check the steel it meant to add is the
    // steel `totalMass()` grew by.
    std::size_t attachedBlocks() const;
    double attachedMass() const;  // kg

    // Global DOF index (3 * node + axis) of each boundary / interior DOF, in the
    // order the reduced model uses. `boundaryDof()` is what a caller couples
    // through and what `staticSolve`'s indices refer to.
    const std::vector<std::uint32_t>& boundaryDof() const;
    const std::vector<std::uint32_t>& interiorDof() const;

    // y = K x over the whole substructure, in global DOF numbering.
    void stiffnessTimes(const std::vector<double>& x, std::vector<double>& y) const;

    // Row `b` of the boundary partition, split: `toBoundary` is row b of K_bb and
    // `toInterior` is row b of K_bi -- which, K being symmetric, is column b of
    // K_ib, and that is what a constraint mode's right-hand side is. Both are
    // resized and cleared. Exposed because forming Psi and the Schur complement
    // through a full sparse multiply per interface DOF costs the whole matrix each
    // time, where the row costs the row.
    void boundaryRow(std::size_t b, std::vector<double>& toBoundary,
                     std::vector<double>& toInterior) const;
    // Row-sum lumped mass per global DOF. Diagonal, so this is the whole of M.
    const std::vector<double>& mass() const;
    double totalMass() const;  // kg; the sum over one axis

    // K_ii^-1 r, in interior numbering, in place. False if the factorisation failed.
    bool interiorSolve(std::vector<double>& r) const;
    // y = K_ii x, in interior numbering.
    void interiorStiffnessTimes(const std::vector<double>& x, std::vector<double>& y) const;

    // The lowest `count` fixed-interface modes, by subspace iteration. Values are
    // eigenvalues of (K_ii, M_ii) -- omega^2 in rad^2/s^2, not frequencies -- and
    // the vectors are **mass-normalised**, phi^T M_ii phi = 1, which is the
    // normalisation that makes the modal mass block the identity. `count` above the
    // interior size is clamped and reported.
    Eigenpairs fixedInterfaceModes(int count, double tolerance = 1e-10,
                                 int maxIterations = 60) const;

    // How many eigenvalues of (K_ii, M_ii) lie strictly below `shift`, by the
    // inertia of an LDL^T factorisation of K_ii - shift * M_ii. Exact -- it counts
    // rather than converges -- and it is what verifies that subspace iteration did
    // not skip a mode. `exact` comes back false if a pivot had to be perturbed,
    // which makes the count a strong estimate rather than a proof.
    int eigenvaluesBelow(double shift, bool* exact = nullptr) const;

    const solidshell::HexMesh& mesh() const;
    const StructuralMaterial& material() const;
    solidshell::Formulation formulation() const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

// --- The reduction --------------------------------------------------------------

struct ReduceParams {
    // Fixed-interface modes to keep. Negative takes the count from
    // `cutoffFrequency` instead, which is the rule §5 recommends.
    int modes = -1;

    // Keep every fixed-interface mode below this, rad/s. 126 rad/s is 20 Hz, twice
    // the 10 Hz top of the hull girder band this engine has loads for -- see §5.
    // The count is settled *exactly*, by an inertia count, before any eigenvector
    // is computed.
    double cutoffFrequency = 126.0;

    // Never keep more than this many however low the cutoff puts them, so a
    // pathological substructure cannot silently cost a dense solve of its own
    // interior. Reported in `Reduction::problems` when it bites.
    int maxModes = 200;

    double eigenTolerance = 1e-10;
    int maxIterations = 60;

    // Run the Sturm sequence check on the mode count. It is one extra banded
    // factorisation; leaving it off is a cost decision, not a correctness one.
    bool verifyModes = true;
};

// The reduced pair, and everything needed to get back to the full model.
//
// Reduced DOF are ordered **boundary first, then modal**: index b in [0, boundary)
// is `Substructure::boundaryDof()[b]`, and index boundary + j is modal coordinate
// j. That ordering is what makes coupling two substructures a matter of matching
// boundary DOF and stacking the modal blocks.
struct Reduction {
    int boundary = 0;
    int modes = 0;
    int size() const { return boundary + modes; }

    // size() x size(), row-major, symmetric to the last bit: the analytically
    // symmetric blocks are computed on one triangle and mirrored.
    std::vector<double> stiffness;  // Pa m -- N/m
    std::vector<double> mass;       // kg

    // Fixed-interface natural frequencies of the modes kept, rad/s, ascending.
    std::vector<double> frequency;
    // The lowest fixed-interface frequency, rad/s, **whether or not it was kept**.
    // Zero modes is only defensible when this is far above the band of interest,
    // so it is reported rather than left to be assumed.
    double firstFixedFrequency = 0;

    // Interior recovery. `constraintModes` is n_i x boundary and `normalModes` is
    // n_i x modes, both row-major over the interior DOF: u_i = Psi u_b + Phi q.
    std::vector<double> constraintModes;
    std::vector<double> normalModes;

    bool modesVerified = false;   // the Sturm check ran and agreed
    int modesBelowCutoff = 0;     // what the inertia count said
    double reduceSeconds = 0;
    std::vector<std::string> problems;

    bool empty() const { return size() == 0; }
};

// Reduce. Always returns something: an unusable substructure yields an empty
// reduction and a full account in `problems`, in the same spirit as
// `makeStructuralMesh` and `buildPatch`.
Reduction craigBampton(const Substructure& substructure, const ReduceParams& params = {});

// --- Using it -------------------------------------------------------------------

// f_r = T^T f, from a load in the substructure's global DOF numbering. Exact for
// any load: the interior part arrives as Psi^T f_i on the boundary and Phi^T f_i
// on the modes.
std::vector<double> reduceLoad(const Substructure& substructure, const Reduction& reduced,
                               const std::vector<double>& load);

// u = T x, in the substructure's global DOF numbering.
std::vector<double> recover(const Substructure& substructure, const Reduction& reduced,
                            const std::vector<double>& state);

// Solve K_r x = f with the reduced DOF listed in `held` fixed at zero. Dense
// Cholesky on the free part. False, with a reason, when what is left is singular --
// which for a free-free substructure with nothing held is exactly right and is what
// the six rigid body modes mean.
bool staticSolve(const Reduction& reduced, const std::vector<double>& load,
                 const std::vector<std::uint32_t>& held, std::vector<double>& state,
                 std::string* problem = nullptr);

// The natural frequencies of the reduced model itself, rad/s ascending, with the
// listed reduced DOF held. This is the quantity property 3 is about: it is an
// upper bound on the full model's, and it falls monotonically as modes are added.
std::vector<double> reducedFrequencies(const Reduction& reduced,
                                       const std::vector<std::uint32_t>& held);

// --- Where the linear model stops ------------------------------------------------

// What the reduced state means in stress, and therefore whether the region is
// still inside the model's validity. See §6: this is a **trigger for promotion**,
// not a strength check, and it under-predicts a concentration.
struct Validity {
    double peakVonMises = 0;   // Pa, over every Gauss point of every element
    double utilisation = 0;    // peakVonMises / yield strength
    int worstElement = -1;
    double peakDisplacement = 0;  // m, largest nodal displacement magnitude
    bool linear = true;           // utilisation < 1: the reduced answer still means something
};

Validity checkValidity(const Substructure& substructure, const Reduction& reduced,
                       const std::vector<double>& state);

// --- Synthesis: coupling two reduced components ----------------------------------
//
// Reduction alone is not component mode *synthesis*. Until two substructures can be
// joined at a shared interface it is one component, and one component is a
// curiosity: the whole point of the tier is that a ship is built from pieces that
// are reduced once and then assembled cheaply.
//
// The coupling is exact and it is nothing more than addition, which is the reason
// the reduced DOF are ordered boundary-first. Two components meeting at an
// interface must have the same displacement there, so their shared boundary DOF
// *are* the same unknown; every other DOF is its own. Assembling is therefore
// scatter-add of both reduced pairs into one, with shared boundary DOF landing on
// the same row and column and the modal blocks stacked side by side. Nothing is
// approximated at this step -- whatever error the assembled model has, it came
// from truncating each component's modes, not from joining them.
//
// **Interface mass is not double counted, and it is worth saying why**, because it
// looks as though it should be. Each element belongs to exactly one component, and
// a component's lumped mass at a shared node comes only from *its own* elements.
// A node on the interface therefore receives part of its mass from each side and
// the sum is the mass the unsplit model had. Splitting a mesh by *node* instead of
// by element would double it, which is the mistake this note exists to prevent.

// Which boundary DOF of two substructures are the same physical DOF.
//
// Matched by position, because two independently built meshes share no numbering.
// A DOF matches only if the node coincides within `tolerance` **and the axis is the
// same** -- without the axis check a coincident node would couple x to y and the
// assembled model would be quietly, plausibly wrong.
struct InterfaceMap {
    // For each boundary index of A, the boundary index of B it coincides with, or
    // -1 for a DOF that A does not share.
    std::vector<int> aToB;
    std::size_t shared = 0;
    double worstGap = 0;  // m, the furthest apart a matched pair actually was
    std::vector<std::string> problems;
};

InterfaceMap matchBoundaries(const Substructure& a, const Substructure& b,
                             double tolerance = 1e-9);

// Two reduced components coupled at their shared boundary.
//
// Assembled DOF are ordered: A's boundary, then B's *unshared* boundary, then A's
// modal coordinates, then B's. `fromA` and `fromB` give the assembled index of
// each component reduced DOF, which is what a caller needs to place a load or read
// a component's state back out.
struct Assembly {
    int boundary = 0;  // the union of the two boundary sets, shared counted once
    int modes = 0;     // a.modes + b.modes
    int size() const { return boundary + modes; }

    std::vector<double> stiffness;  // N/m
    std::vector<double> mass;       // kg

    std::vector<int> fromA;  // size = reduction A's size()
    std::vector<int> fromB;  // size = reduction B's size()

    std::vector<std::string> problems;
    bool empty() const { return size() == 0; }
};

// **Two components, and only two.** There is no way to join an `Assembly` to a
// third component, because `InterfaceMap` is expressed in the *substructures'*
// boundary DOF and an `Assembly` has no substructure behind it. A ship is many
// components, so this is a real limit and not a simplification for exposition.
//
// Generalising it wants one change rather than a rewrite: an assembled model needs
// to carry the boundary DOF identity its components had -- which global DOF of
// which mesh each assembled boundary row *is* -- so that the same position match
// can be run against it. The scatter-add underneath already does not care how many
// components it is given.
//
// **There is now a pair worth assembling**, which there was not when the note
// above was written: `coupling.{hpp,cpp}` joins a Tier-2 zone to the plating round
// it and drives the zone's perimeter from the result. Two components is exactly
// enough for that, so the three-component limit is no longer the binding one.
//
// **And the mesher exists now too** -- `section.{hpp,cpp}` cuts a region between two
// transverse planes and hands it over with `nodesNearPlanes`' interface and an
// `Attachment` of its longitudinals. What it cannot do is weld a junction, so a
// hold of the reference ferry arrives as seven disconnected pieces; that is a
// property of `makeStructuralMesh`, which shares no corner between two panel roles,
// and of the solid-shell, whose node pair carries one thickness direction where a
// corner has two. The consequence for this file is worth knowing before reducing
// one: the section's lowest fixed-interface frequency is its *decks*' -- 0.78 Hz on
// that hold, against 3.46 Hz for the shell alone -- so the 20 Hz `cutoffFrequency`
// below asks for 178 modes, takes 275 s, and does not converge. Guyan is 6 s and is
// exactly right at the interface.
Assembly assemble(const Reduction& a, const Reduction& b, const InterfaceMap& map);

// Natural frequencies of an assembled model, rad/s ascending, with the listed
// assembled DOF held. The same upper-bound property holds as for one component.
std::vector<double> assembledFrequencies(const Assembly& assembly,
                                         const std::vector<std::uint32_t>& held);

// Solve K x = f on an assembled model, with `held` fixed at zero.
bool assembledStaticSolve(const Assembly& assembly, const std::vector<double>& load,
                          const std::vector<std::uint32_t>& held, std::vector<double>& state,
                          std::string* problem = nullptr);

// The part of an assembled state belonging to one component, in that component's
// own reduced DOF order -- what `recover()` wants. Empty if `state` is not this
// assembly's: every index is otherwise in range, so a short state would come back
// as a plausible field quietly missing its modal content.
std::vector<double> componentState(const Assembly& assembly, const std::vector<int>& from,
                                   const std::vector<double>& state);

}  // namespace sim::reduction
