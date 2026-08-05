// SPDX-License-Identifier: MIT
//
// Validation of the multi-point constraint and the eccentric stiffener built on
// it -- `engine/sim/constraint.{hpp,cpp}`.
//
// The claim being tested is `zone.hpp` §3's: that a solid-shell's two nodes
// through the thickness *are* the rotation of the plate's cross-section, so an
// eccentric member is an exact linear function of that pair and needs no
// rotational degree of freedom anywhere. Three things have to be true for that to
// be worth anything, and each has its own closed form:
//
//   * the tie is **exact**, including for a finite rigid rotation, or a stiffener
//     on a patch that has swung would grow force out of nothing;
//   * the eccentricity is **right**, which means the stiffened panel's second
//     moment about the *combined* neutral axis has to come out of the finite
//     element model equal to what `scantlings::stiffenedSection` computes from the
//     section alone -- two independent routes to one number, and the number this
//     whole exercise exists for;
//   * the tripping mode is **not zero energy**, which is the specific defect the
//     hinge formulation has and which no bending test can see.
//
// The section-property checks are run through the *linear* assembly rather than
// through `zone::Solver`, and that is deliberate. The pure-bending field is an
// exact solution of elasticity, so a linear assembly reproduces its energy to
// rounding and the comparison is an identity; running it through the explicit
// solver adds the co-rotational formulation's own O(theta^2) error, which at
// kappa = 1e-3 over a two-metre patch is 1.6% and would turn an identity into a
// tolerance. `test_zone.cpp` makes the same comparison through the solver, where
// that 1.6% is the thing being measured.
#include "engine/sim/constraint.hpp"
#include "engine/sim/plasticity.hpp"
#include "engine/sim/reduction.hpp"
#include "engine/sim/scantlings.hpp"
#include "engine/sim/solid_shell.hpp"
#include "engine/sim/zone.hpp"
#include "harness.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <string>
#include <vector>

using namespace sim;
using testing::expectEqual;
using testing::expectNear;
using testing::expectTrue;

namespace {

void expectEqualCount(const std::string& what, std::size_t got, std::size_t want) {
    testing::expectEqual(what, static_cast<long long>(got), static_cast<long long>(want));
}

// --- Fixtures ------------------------------------------------------------------

// A flat rectangle of plating, `nx` by `ny` panels, optionally with one
// longitudinal down the middle. The member is on the panel seam at y = 0, which is
// where `makeStructuralMesh` puts them, so the mesher's grid points land on it
// exactly.
StructuralMesh flatStrip(double lengthX, double spanY, double thickness, int nx, int ny,
                         const StiffenerProfile& profile, bool stiffened) {
    StructuralMesh mesh;
    mesh.materials.push_back(ah36Steel());
    mesh.frameSpacing = lengthX;
    for (int i = 0; i < nx; ++i)
        for (int j = 0; j < ny; ++j) {
            PlatePanel p;
            const double x0 = -0.5 * lengthX + lengthX * i / nx;
            const double x1 = -0.5 * lengthX + lengthX * (i + 1) / nx;
            const double y0 = -0.5 * spanY + spanY * j / ny;
            const double y1 = -0.5 * spanY + spanY * (j + 1) / ny;
            p.corner[0] = {x0, y0, 0};
            p.corner[1] = {x1, y0, 0};
            p.corner[2] = {x1, y1, 0};
            p.corner[3] = {x0, y1, 0};
            p.thickness = thickness;
            p.material = 0;
            p.role = PanelRole::Shell;
            mesh.panels.push_back(p);
        }
    if (stiffened) {
        StructuralMember member;
        member.a = {-0.5 * lengthX, 0.0, 0.0};
        member.b = {0.5 * lengthX, 0.0, 0.0};
        // The patch's outward normal is +z, and a shell longitudinal's web rises
        // *inward* -- `Section::at` hands out the inward normal for exactly this.
        // So the eccentricity is negative along the node pair, and every sign in
        // the tie has to survive that.
        member.rise = {0, 0, -1};
        member.profile = profile;
        member.attachedPlateThickness = thickness;
        member.role = MemberRole::Longitudinal;
        mesh.members.push_back(member);
    }
    return mesh;
}

zone::Patch stiffenedPatch(double lengthX, double spanY, double thickness, int nx, int ny,
                           const StiffenerProfile& profile, bool stiffened, int subdivision,
                           zone::Edge edge = zone::Edge::Free) {
    const StructuralMesh strip = flatStrip(lengthX, spanY, thickness, nx, ny, profile, stiffened);
    zone::MeshParams params;
    params.radius = 1000.0;  // the whole strip, so the geometry is the one asked for
    params.subdivision = subdivision;
    params.stiffeners = stiffened ? zone::Stiffeners::Modelled : zone::Stiffeners::Ignored;
    params.outward = {0, 0, 1};
    params.edge = edge;
    return zone::buildPatch(strip, {0, 0, 0}, params);
}

// The exact pure-bending elasticity field about `axis`, evaluated at every node.
// sigma_xx = E kappa (z - axis) and nothing else, which is the state
// `zone::Preload` imposes and the one whose stored energy is E kappa^2 I / 2 per
// unit length about that axis, for *any* axis, by the parallel axis theorem.
std::vector<double> pureBending(const solidshell::HexMesh& mesh, double kappa, double axis,
                                double poisson) {
    std::vector<double> field(mesh.nodeCount() * 3, 0.0);
    for (std::size_t node = 0; node < mesh.nodeCount(); ++node) {
        const double x = mesh.position[node * 3];
        const double y = mesh.position[node * 3 + 1];
        const double h = mesh.position[node * 3 + 2] - axis;
        field[node * 3] = kappa * h * x;
        field[node * 3 + 1] = -poisson * kappa * h * y;
        field[node * 3 + 2] =
            -poisson * 0.5 * kappa * h * h - 0.5 * kappa * x * x + 0.5 * poisson * kappa * y * y;
    }
    return field;
}

// Strain energy of a nodal field, plating and fibres separately: u^T K u / 2 with
// the element stiffness and the fibres' condensed rank-one blocks.
struct FieldEnergy { double plating = 0, fibre = 0; double total() const { return plating + fibre; } };

FieldEnergy fieldEnergy(const zone::Patch& patch, const constraint::Stiffening& stiffening,
                        const std::vector<double>& field) {
    FieldEnergy out;
    for (std::size_t e = 0; e < patch.elementCount(); ++e) {
        double nodes[solidshell::kDof], u[solidshell::kDof];
        patch.mesh.gather(e, patch.mesh.position, nodes);
        patch.mesh.gather(e, field, u);
        double k[solidshell::kDof * solidshell::kDof];
        solidshell::elementStiffness(nodes, patch.material, solidshell::Formulation::SolidShell, k);
        for (int i = 0; i < solidshell::kDof; ++i) {
            double sum = 0;
            for (int j = 0; j < solidshell::kDof; ++j) sum += k[i * solidshell::kDof + j] * u[j];
            out.plating += 0.5 * u[i] * sum;
        }
    }
    const constraint::RestFibers forms = constraint::restFibers(stiffening, patch.mesh.position);
    for (std::size_t f = 0; f < stiffening.fiberCount(); ++f) {
        const constraint::FiberStiffness block = constraint::fiberStiffness(
            stiffening.fiber[f], patch.mesh.position, forms.length[f],
            patch.material.youngsModulus);
        if (!(block.scale > 0)) continue;
        double projection = 0;
        for (int i = 0; i < 12; ++i) projection += block.vector[i] * field[block.dof[i]];
        out.fibre += 0.5 * block.scale * projection * projection;
    }
    return out;
}

// The whole stiffness, densely, for an eigenvalue question. Small meshes only.
std::vector<double> denseStiffness(const zone::Patch& patch,
                                   const constraint::Stiffening& stiffening) {
    const std::size_t n = patch.nodeCount() * 3;
    std::vector<double> k(n * n, 0.0);
    for (std::size_t e = 0; e < patch.elementCount(); ++e) {
        double nodes[solidshell::kDof];
        patch.mesh.gather(e, patch.mesh.position, nodes);
        double ke[solidshell::kDof * solidshell::kDof];
        solidshell::elementStiffness(nodes, patch.material, solidshell::Formulation::SolidShell, ke);
        std::size_t dof[solidshell::kDof];
        for (int a = 0; a < solidshell::kNodes; ++a) {
            const std::size_t node =
                patch.mesh.index[e * solidshell::kNodes + static_cast<std::size_t>(a)];
            for (int i = 0; i < 3; ++i) dof[a * 3 + i] = node * 3 + static_cast<std::size_t>(i);
        }
        for (int p = 0; p < solidshell::kDof; ++p)
            for (int q = 0; q < solidshell::kDof; ++q)
                k[dof[p] * n + dof[q]] += ke[p * solidshell::kDof + q];
    }
    const constraint::RestFibers forms = constraint::restFibers(stiffening, patch.mesh.position);
    for (std::size_t f = 0; f < stiffening.fiberCount(); ++f) {
        const constraint::FiberStiffness block = constraint::fiberStiffness(
            stiffening.fiber[f], patch.mesh.position, forms.length[f],
            patch.material.youngsModulus);
        if (!(block.scale > 0)) continue;
        for (int p = 0; p < 12; ++p)
            for (int q = 0; q < 12; ++q)
                k[static_cast<std::size_t>(block.dof[p]) * n + block.dof[q]] +=
                    block.scale * block.vector[p] * block.vector[q];
    }
    return k;
}

// --- 1. The tie -----------------------------------------------------------------
//
// The whole formulation rests on one claim: a point at through-thickness offset
// `e` is an exact linear function of the shell's node pair, with no rotational
// degree of freedom anywhere. Two identities say so, and they are identities and
// not tolerances.

void testTheTieIsExactUnderRigidMotionAndUnderRotationOfThePlate() {
    std::printf("\n   the tie: a linear function of the pair, exact for a finite rotation\n");
    const double thickness = 0.012;

    // The weights, against their own closed form. The bottom face, the top face,
    // the mid-surface, and the extrapolation an eccentric stiffener needs.
    expectNear("the -zeta face has weight zero", constraint::tieWeight(-0.5 * thickness, thickness),
               0.0, 1e-15);
    expectNear("the +zeta face has weight one", constraint::tieWeight(0.5 * thickness, thickness),
               1.0, 1e-15);
    expectNear("the mid-surface has weight one half", constraint::tieWeight(0.0, thickness), 0.5,
               1e-15);
    const double far = constraint::tieWeight(0.1, thickness);
    expectNear("and 100 mm off 12 mm plating extrapolates", far, (0.1 + 0.006) / 0.012, 1e-14);
    expectTrue("which is an extrapolation, not an interpolation", far > 8.0);

    // One pair, and a tied point 100 mm off the mid-surface.
    std::vector<double> nodal(6, 0.0);
    const Vec3 bottom{0.3, -0.2, -0.5 * thickness}, top{0.3, -0.2, 0.5 * thickness};
    for (int k = 0; k < 3; ++k) {
        nodal[static_cast<std::size_t>(k)] = bottom[k];
        nodal[3 + static_cast<std::size_t>(k)] = top[k];
    }
    constraint::Tie tie{0, 1, constraint::tieWeight(-0.1, thickness)};
    const Vec3 rest = constraint::tiedPoint(tie, nodal);
    expectNear("the tied point sits at the offset it was given", rest.z, -0.1, 1e-15);
    expectNear("directly above its pair in x", rest.x, bottom.x, 1e-15);
    expectNear("and in y", rest.y, bottom.y, 1e-15);

    // **A finite rigid rotation.** 0.7 rad about an oblique axis, plus a
    // translation. The tie is linear, so the tied point has to land exactly on
    // R x + c -- if it did not, a stiffener on a swinging patch would strain.
    const double angle = 0.7;
    const Vec3 axis = normalize(Vec3{0.3, -0.6, 0.74});
    const Vec3 shift{1.5, -2.5, 0.75};
    const auto rotate = [&](const Vec3& p) {
        const double c = std::cos(angle), s = std::sin(angle);
        return p * c + cross(axis, p) * s + axis * (dot(axis, p) * (1.0 - c));
    };
    std::vector<double> moved(6, 0.0);
    for (int end = 0; end < 2; ++end) {
        const Vec3 turned = rotate(end == 0 ? bottom : top) + shift;
        for (int k = 0; k < 3; ++k)
            moved[static_cast<std::size_t>(end) * 3 + static_cast<std::size_t>(k)] = turned[k];
    }
    const Vec3 got = constraint::tiedPoint(tie, moved);
    const Vec3 want = rotate(rest) + shift;
    const double error = length(got - want);
    std::printf("     0.70 rad about (%.2f, %.2f, %.2f): the tied point lands %.3e m from"
                " R x + c, having travelled %.4f m\n", axis.x, axis.y, axis.z, error,
                length(want - rest));
    // 1e-14 rather than 1e-9, and the bound is derived rather than picked. The tie
    // is `9.33 x_bottom - 8.33 x_top` at this offset, so it amplifies the last bit
    // of a metre-sized coordinate by `1 + 2|w| = 18`; 18 ulp of 3 m is 1.2e-14.
    // Measured at 3.6e-15 -- under an ulp of the amplified operand, which is the
    // statement that the tie *is* the linear map and not merely close to one.
    expectTrue("a finite rigid rotation carries the tied point exactly", error < 1e-14);
    // Non-vacuous: the rotation actually moved it, so "exact" is not "nothing
    // happened".
    expectTrue("and the rotation was a real one", length(want - rest) > 0.5);

    // **The plate rotating about its own mid-surface**, which is the motion an
    // eccentric member exists to feel. Separating the pair by t*phi across the
    // seam tilts the through-thickness fibre by phi, and the tied point at offset
    // e has to move by e*phi. This is the eccentric kinematics, and it is what a
    // model that tied the member to the mid-surface would get identically zero.
    // Applied to the *displacement* rather than to the position, so the answer is
    // not the difference of two 0.2 m numbers and the assertion can be an identity.
    const double phi = 1.0e-3;
    std::vector<double> tip(6, 0.0);
    tip[1] = -0.5 * thickness * phi;
    tip[4] = 0.5 * thickness * phi;
    expectNear("a plate rotation of phi moves a member at offset e by e phi",
               constraint::tiedPoint(tie, tip).y, -0.1 * phi, 1e-19);
    const constraint::Tie flat{0, 1, 0.5};
    expectTrue("and a member tied to the mid-surface does not move at all",
               constraint::tiedPoint(flat, tip).y == 0.0);

    // The force distribution is the transpose of the interpolation, which is what
    // makes the condensed system symmetric. Virtual work says so: f . u(tied) has
    // to equal the distributed force dotted with the master displacements, for any
    // displacement at all.
    const Vec3 force{3.0e5, -1.2e6, 7.5e5};
    std::vector<double> probe{0.11, -0.23, 0.37, -0.41, 0.53, 0.19};
    std::vector<double> distributed(6, 0.0);
    constraint::applyTiedForce(tie, force, distributed);
    double work = 0;
    for (std::size_t d = 0; d < 6; ++d) work += distributed[d] * probe[d];
    const Vec3 at = constraint::tiedPoint(tie, probe);
    expectNear("the force distribution is the transpose of the interpolation", work, dot(force, at),
               1e-12 * std::abs(work));
    expectTrue("and the work was not zero", std::abs(work) > 1.0e4);
}

// --- 2. The fibres carry the profile ---------------------------------------------

void testTheFibresReproduceTheProfileSection() {
    std::printf("\n   profile fibres against scantlings::profileSection\n");
    const double thickness = 0.012;
    const StiffenerProfile shapes[] = {flatBar(0.200, 0.010), flatBar(0.120, 0.008),
                                       tee(0.300, 0.011, 0.150, 0.016),
                                       angle(0.250, 0.010, 0.090, 0.014)};
    const char* names[] = {"200x10 flat bar", "120x8 flat bar", "300x11+150x16 tee",
                           "250x10+90x14 angle"};

    for (int shape = 0; shape < 4; ++shape) {
        const ProfileSection section = profileSection(shapes[shape]);
        // Both senses, because a shell longitudinal rises against the node pair and
        // a deck one may rise with it, and a sign error here is a wrong neutral
        // axis that every symmetric test would miss.
        for (double sign : {-1.0, 1.0}) {
            const constraint::ProfileFibers fibres =
                constraint::profileFibers(shapes[shape], thickness, sign);
            double area = 0, first = 0, second = 0;
            for (int i = 0; i < fibres.count; ++i) {
                // Back to the height above the attachment face, which is what
                // `profileSection` measures from.
                const double height = sign * fibres.offset[i] - 0.5 * thickness;
                area += fibres.area[i];
                first += fibres.area[i] * height;
                second += fibres.area[i] * height * height;
            }
            const double aboutFace =
                section.secondMoment + section.area * section.centroid * section.centroid;
            if (sign < 0)
                std::printf("     %-20s %d fibres: A %.6e  c %.6f  I_face %.6e  (I_own is %.0f%%"
                            " of it, and is what the second Gauss station buys)\n", names[shape],
                            fibres.count, area, first / area, second,
                            100.0 * section.secondMoment /
                                (section.secondMoment +
                                 section.area * section.centroid * section.centroid));
            expectNear(std::string(names[shape]) + ": the fibres carry the profile's area", area,
                       section.area, 1e-15 * section.area);
            expectNear(std::string(names[shape]) + ": and its centroid", first / area,
                       section.centroid, 1e-14 * section.centroid);
            // The second moment is the one two-point Gauss buys and one-point does
            // not, and it is the difference between carrying `I_own` and losing it.
            expectNear(std::string(names[shape]) + ": and its second moment about the face",
                       second, aboutFace, 1e-14 * aboutFace);
        }
        // Non-vacuous: the profile's *own* second moment has to be a real share of
        // the total, or the identity above would be satisfied by a one-station rule
        // that only got the area and the first moment right.
        const double aboutFace =
            section.secondMoment + section.area * section.centroid * section.centroid;
        expectTrue(std::string(names[shape]) + ": I_own is a real share of I about the face",
                   section.secondMoment > 0.15 * aboutFace);
    }
}

// --- 3. The stiffened panel, against its own section modulus ---------------------
//
// The measurement this whole exercise exists for. `scantlings.hpp` §1 rejects
// smearing because the panel's second moment falls by 130x at identical area, so
// the eccentricity is the point and a formulation that gets the offset wrong will
// look plausible and be badly wrong.
//
// The instrument is the parallel axis theorem. Impose pure bending about an axis
// `z0` and the stored energy is `E kappa^2 (I_NA + A (z0 - z_NA)^2) L / 2` -- a
// parabola in z0 whose *minimum* is the neutral axis, whose *value there* is the
// second moment, and whose *curvature* is the area. Sweeping z0 therefore checks
// all three of `stiffenedSection`'s outputs against the finite element model, and
// an offset that is wrong moves the vertex.

void testTheStiffenedPanelMatchesItsSectionProperties() {
    std::printf("\n   the stiffened panel's second moment, against stiffenedSection\n");
    const double lengthX = 2.0, spanY = 0.7, thickness = 0.012, kappa = 1.0e-3;
    const StiffenerProfile bar = flatBar(0.200, 0.010);
    const StructuralMaterial steel = ah36Steel();
    const StiffenedSection section = stiffenedSection(bar, thickness, spanY);
    // `stiffenedSection` measures towards the stiffener; the patch's +z is away
    // from it, so the neutral axis is at -neutralAxis in the patch's own frame.
    const double neutralAxis = -section.neutralAxis;

    const zone::Patch bare = stiffenedPatch(lengthX, spanY, thickness, 4, 2, bar, false, 2);
    const zone::Patch stiff = stiffenedPatch(lengthX, spanY, thickness, 4, 2, bar, true, 2);
    expectTrue("both patches meshed", !bare.empty() && !stiff.empty());
    expectTrue("the stiffened one carries fibres", stiff.stiffening.fiberCount() > 0);
    expectEqual("from one member", static_cast<long long>(stiff.stiffening.members), 1LL);
    expectNear("covering the whole seam", stiff.stiffening.length, lengthX, 1e-12 * lengthX);
    expectNear("with the profile's own mass", stiff.stiffening.mass,
               profileSection(bar).area * lengthX * steel.density,
               1e-12 * profileSection(bar).area * lengthX * steel.density);

    std::printf("     %-10s %14s %14s %14s %12s\n", "axis (m)", "plating (J)", "fibres (J)",
                "closed form", "ratio");
    for (double axis : {0.0, neutralAxis, -0.05, 0.05, 0.15}) {
        const std::vector<double> field =
            pureBending(stiff.mesh, kappa, axis, steel.poissonRatio);
        const FieldEnergy got = fieldEnergy(stiff, stiff.stiffening, field);
        const double offset = axis - neutralAxis;
        const double second = section.secondMoment + section.area * offset * offset;
        const double want = 0.5 * steel.youngsModulus * kappa * kappa * second * lengthX;
        std::printf("     %-10.6f %14.6e %14.6e %14.6e %12.9f\n", axis, got.plating, got.fibre,
                    want, got.total() / want);
        // 1e-8, and measured at 2e-10. Both routes are exact: the pure-bending
        // field is an elasticity solution the element reproduces to rounding (a
        // single element does, to ten digits), and two-point Gauss through each
        // rectangle of the profile integrates the fibres' quadratic energy density
        // exactly. Anything looser here would pass on a model that had the area
        // right and the eccentricity approximately right, which is exactly the
        // failure being guarded against.
        expectNear("the stiffened panel's stored energy is E kappa^2 I L / 2", got.total(), want,
                   1e-8 * want);
    }

    // The vertex of that parabola *is* the neutral axis, recovered from three
    // measurements rather than read off the model. This is the check an offset
    // error moves: a stiffener attached at the mid-surface puts it at zero.
    const auto energyAt = [&](double axis) {
        return fieldEnergy(stiff, stiff.stiffening,
                           pureBending(stiff.mesh, kappa, axis, steel.poissonRatio))
            .total();
    };
    const double step = 0.02;
    const double left = energyAt(-step), centre = energyAt(0.0), right = energyAt(step);
    const double vertex = -0.5 * step * (right - left) / (right - 2.0 * centre + left);
    const double curvature = (right - 2.0 * centre + left) / (step * step);
    const double areaFromCurvature = curvature / (steel.youngsModulus * kappa * kappa * lengthX);
    std::printf("     the parabola through three points: vertex %.8f m (want %.8f), area from its"
                " curvature %.8e m^2 (want %.8e)\n", vertex, neutralAxis, areaFromCurvature,
                section.area);
    expectNear("the neutral axis is where the parabola turns", vertex, neutralAxis,
               1e-9 * std::abs(neutralAxis));
    expectNear("and its curvature is the section's area", areaFromCurvature, section.area,
               1e-9 * section.area);

    // --- The negative controls ----------------------------------------------
    //
    // Without them the checks above pass on a constraint that does nothing.
    const double bareEnergy =
        fieldEnergy(bare, bare.stiffening, pureBending(bare.mesh, kappa, 0.0, steel.poissonRatio))
            .total();
    const double stiffEnergy = energyAt(0.0);
    const double plateSecond = spanY * thickness * thickness * thickness / 12.0;
    expectNear("the bare plate is its own b t^3 / 12", bareEnergy,
               0.5 * steel.youngsModulus * kappa * kappa * plateSecond * lengthX,
               1e-8 * 0.5 * steel.youngsModulus * kappa * kappa * plateSecond * lengthX);
    const double gain = stiffEnergy / bareEnergy;
    std::printf("     the stiffener multiplies the panel's bending stiffness by %.1f about the"
                " plate's own mid-surface, and by %.1f about each section's own neutral axis\n",
                gain, section.secondMoment / plateSecond);
    expectTrue("the stiffener is not a rounding effect", gain > 100.0);

    // **Eccentricity zeroed**: the same fibres, the same area, tied to the
    // mid-surface instead. That is the model a wrong tie weight produces, and it is
    // `scantlings.hpp` §1's smearing argument written as an energy: with every
    // fibre at one point the *whole* stiffener drops out of the bending stiffness,
    // Steiner term and own second moment together, and the panel is exactly the
    // bare plate again. Not "softer" -- identical, which is a much stronger thing
    // to be able to assert.
    constraint::Stiffening flat = stiff.stiffening;
    for (constraint::Fiber& fibre : flat.fiber)
        for (constraint::Tie& tie : fibre.end) tie.weight = 0.5;
    const double flatEnergy =
        fieldEnergy(stiff, flat, pureBending(stiff.mesh, kappa, 0.0, steel.poissonRatio)).total();
    std::printf("     the same fibres with no eccentricity: %.6e J against the bare plate's"
                " %.6e and the eccentric tie's %.6e\n", flatEnergy, bareEnergy, stiffEnergy);
    expectNear("tying the fibres to the mid-surface leaves the bare plate and nothing else",
               flatEnergy, bareEnergy, 1e-12 * bareEnergy);

    // **One fibre per rectangle instead of two**: the same area at the same
    // centroid, which gets the Steiner term right and loses the profile's own
    // second moment. That is the arithmetic error two-point Gauss exists to
    // prevent, and for this bar `I_own` is 27% of the stiffened section, so the
    // model would come out a quarter soft while looking entirely plausible.
    constraint::Stiffening lumped;
    lumped.material = stiff.stiffening.material;
    const double centroidOffset = -(0.5 * thickness + profileSection(bar).centroid);
    for (std::size_t f = 0; f < stiff.stiffening.fiberCount(); f += 2) {
        constraint::Fiber merged = stiff.stiffening.fiber[f];
        merged.area = stiff.stiffening.fiber[f].area + stiff.stiffening.fiber[f + 1].area;
        merged.offset = centroidOffset;
        for (constraint::Tie& tie : merged.end)
            tie.weight = constraint::tieWeight(centroidOffset, thickness);
        lumped.fiber.push_back(merged);
    }
    const double lumpedEnergy = fieldEnergy(stiff, lumped,
                                            pureBending(stiff.mesh, kappa, 0.0, steel.poissonRatio))
                                    .total();
    const double lumpedSecond =
        section.secondMoment + section.area * neutralAxis * neutralAxis -
        profileSection(bar).secondMoment;
    std::printf("     one fibre per rectangle instead of two: %.6e J against %.6e, %.1f%% soft,"
                " and losing I_own predicts %.1f%%\n", lumpedEnergy, stiffEnergy,
                100.0 * (1.0 - lumpedEnergy / stiffEnergy),
                100.0 * profileSection(bar).secondMoment /
                    (section.secondMoment + section.area * neutralAxis * neutralAxis));
    expectNear("a single station loses exactly the profile's own second moment", lumpedEnergy,
               0.5 * steel.youngsModulus * kappa * kappa * lumpedSecond * lengthX,
               1e-8 * 0.5 * steel.youngsModulus * kappa * kappa * lumpedSecond * lengthX);
    expectTrue("which is a fifth of the panel, not a rounding difference",
               lumpedEnergy < 0.85 * stiffEnergy);
}

// --- 4. The tripping mode --------------------------------------------------------
//
// `zone.hpp` §3's first objection to a meshed web is that it is a hinge: it keeps
// an independent rotation about the seam, that rotation costs no energy, and an
// explicit scheme turns a zero-energy mode into ringing rather than into an
// answer. The fix here is structural rather than numerical -- the member has no
// degrees of freedom of its own at all -- so the statement to check is that the
// stiffened system's null space is exactly the plating's.

void testTrippingIsNotAZeroEnergyMode() {
    std::printf("\n   the tripping mode: the stiffened system's null space, and what it costs\n");
    const double side = 0.7, thickness = 0.012, phi = 1.0e-3;
    const StiffenerProfile bar = flatBar(0.200, 0.010);

    const zone::Patch bare = stiffenedPatch(side, side, thickness, 2, 2, bar, false, 2);
    const zone::Patch stiff = stiffenedPatch(side, side, thickness, 2, 2, bar, true, 2);
    const int dof = static_cast<int>(stiff.nodeCount() * 3);
    expectEqual("the two patches have the same degrees of freedom",
                static_cast<long long>(bare.nodeCount()), static_cast<long long>(stiff.nodeCount()));
    expectTrue("and the stiffened one has fibres", stiff.stiffening.fiberCount() > 0);

    const auto zeroModes = [&](const zone::Patch& patch, const constraint::Stiffening& stiffening,
                               const char* label) {
        const reduction::Eigenpairs pairs =
            reduction::symmetricEigen(denseStiffness(patch, stiffening), dof);
        expectTrue(std::string(label) + ": the spectrum converged", pairs.converged);
        // **The threshold is derived, not chosen.** This repo has twice been bitten
        // by a fixed rigid-body cutoff landing between the translations and the
        // rotations; here the elastic scale is the first non-zero eigenvalue, and
        // anything a millionth of it is numerical zero.
        double elastic = 0;
        for (int i = 0; i < dof; ++i)
            if (pairs.value[static_cast<std::size_t>(i)] > 1.0) {
                elastic = pairs.value[static_cast<std::size_t>(i)];
                break;
            }
        int zeros = 0;
        for (int i = 0; i < dof; ++i)
            if (std::abs(pairs.value[static_cast<std::size_t>(i)]) < 1e-6 * elastic) ++zeros;
        std::printf("     %-10s first elastic eigenvalue %.4e N/m, largest |lambda| below it"
                    " %.3e, zero modes %d\n", label, elastic,
                    std::abs(pairs.value[static_cast<std::size_t>(zeros - 1)]), zeros);
        return zeros;
    };
    const int bareZeros = zeroModes(bare, bare.stiffening, "plating");
    const int stiffZeros = zeroModes(stiff, stiff.stiffening, "stiffened");
    expectEqual("a free plating patch has exactly the six rigid body modes", bareZeros, 6LL);
    // The whole claim, in one line: a member condensed onto the plating's own
    // degrees of freedom cannot add a mechanism, because it adds no degrees of
    // freedom. The hinge is exactly the formulation that does.
    expectEqual("and the stiffened patch has exactly the same six, and no seventh", stiffZeros,
                6LL);

    // The tripping motion itself: every pair on the seam separated across it, which
    // tips the web about the seam by phi.
    std::vector<double> tip(stiff.nodeCount() * 3, 0.0);
    int stations = 0;
    for (std::size_t node = 0; node + 1 < stiff.nodeCount(); node += 2) {
        if (std::abs(stiff.mesh.position[node * 3 + 1]) > 1e-9) continue;
        ++stations;
        tip[node * 3 + 1] = -0.5 * thickness * phi;
        tip[(node + 1) * 3 + 1] = 0.5 * thickness * phi;
    }
    expectTrue("the seam has stations to tip", stations >= 3);

    const FieldEnergy tipped = fieldEnergy(stiff, stiff.stiffening, tip);
    const FieldEnergy tippedBare = fieldEnergy(stiff, constraint::Stiffening{}, tip);
    std::printf("     tipping the web by %.0e rad costs %.6e J; the plating alone costs %.6e J\n",
                phi, tipped.total(), tippedBare.total());
    expectTrue("the tripping motion costs energy", tipped.total() > 0);
    // Non-vacuous: "positive" has to mean positive on the scale of the structure,
    // not on the scale of rounding. The plating's own bending energy under this
    // motion is the scale.
    expectTrue("and it is a structural stiffness, not a rounding residue",
               tipped.total() > 1.0e-3 * phi * phi * 206.0e9 * thickness * thickness * thickness);

    // **And the honest half of the answer.** A bar has stiffness along its own axis
    // and none across it, so the fibres contribute *exactly* nothing to a motion
    // that only tips them sideways. What restrains tripping here is the plating,
    // through the tie; the web is forced to follow the plate's cross-section, which
    // is the opposite error from the hinge -- over-restrained rather than free.
    // Asserted as an identity because "small" would leave it open whether the
    // fibres were contributing a little and wrongly.
    expectTrue("the fibres contribute exactly nothing to tripping",
               tipped.fibre == 0.0 && tipped.total() == tippedBare.total());

    // The tie carries the eccentric kinematics under that motion, exactly: the
    // outermost fibre moves by its own offset times the rotation. This is what
    // "the constraint is enforced" means here -- the tied point has no degrees of
    // freedom to drift, so the check has to be the physical consequence rather
    // than the algebraic identity, which would be vacuous.
    double worst = 0, largest = 0;
    for (const constraint::Fiber& fibre : stiff.stiffening.fiber)
        for (const constraint::Tie& tie : fibre.end) {
            const double moved = constraint::tiedPoint(tie, tip).y;
            worst = std::max(worst, std::abs(moved - fibre.offset * phi));
            largest = std::max(largest, std::abs(moved));
        }
    std::printf("     every fibre moves laterally by its own offset times phi, to %.3e m of"
                " a largest motion of %.3e m\n", worst, largest);
    expectTrue("the tie carries the eccentricity to rounding", worst < 1e-18);
    expectTrue("and the motion it carries is not zero", largest > 1e-5);

    // Contrast, so "the fibres do nothing" is a statement about *this* motion and
    // not about the fibres: under bending they do almost all of it.
    const FieldEnergy bent =
        fieldEnergy(stiff, stiff.stiffening, pureBending(stiff.mesh, 1e-3, 0.0, 0.3));
    expectTrue("the same fibres carry the bending mode", bent.fibre > 10.0 * bent.plating);
}

// The number the tripping restraint actually is, against a closed form: a plate
// strip of width b, seam down the middle, clamped at both far edges, given a
// rotation phi at the seam, resists at `16 D / b` per unit length -- two strips of
// width b/2, each `4 D / (b/2)`. Cylindrical bending, so the modulus is the plate
// rigidity and not E I, which is what pinning the length direction imposes.
void testTheSeamRotationStiffnessIsThePlatings() {
    std::printf("\n   the seam's rotational stiffness against 16 D / b\n");
    const StructuralMaterial steel = ah36Steel();
    const double thickness = 0.012, span = 0.7, run = 0.35, phi = 1.0e-3;
    const double rigidity = steel.youngsModulus * thickness * thickness * thickness /
                            (12.0 * (1.0 - steel.poissonRatio * steel.poissonRatio));
    const double want = 0.5 * (16.0 * rigidity / span) * phi * phi * run;

    double last = 0;
    for (int refine : {1, 2, 4}) {
        solidshell::HexMesh mesh =
            solidshell::makePlateMesh(run, span, thickness, refine, 8 * refine, 1);
        for (std::size_t node = 0; node < mesh.nodeCount(); ++node) mesh.pin(node, 0, 0.0);
        for (std::size_t node = 0; node < mesh.nodeCount(); ++node) {
            const double y = mesh.position[node * 3 + 1];
            const double z = mesh.position[node * 3 + 2];
            if (std::abs(y) < 1e-9 || std::abs(y - span) < 1e-9) {
                mesh.pin(node, 1, 0.0);
                mesh.pin(node, 2, 0.0);
            } else if (std::abs(y - 0.5 * span) < 1e-9) {
                mesh.pin(node, 1, (z > 0 ? 0.5 : -0.5) * thickness * phi);
                mesh.pin(node, 2, 0.0);
            }
        }
        std::vector<double> load(mesh.nodeCount() * 3, 0.0), displacement;
        std::string problem;
        expectTrue("the seam rotation solves: " + problem,
                   solidshell::solveStatic(mesh, steel, solidshell::Formulation::SolidShell, load,
                                           displacement, &problem));
        double energy = 0;
        for (std::size_t e = 0; e < mesh.elementCount(); ++e) {
            double nodes[solidshell::kDof], u[solidshell::kDof];
            mesh.gather(e, mesh.position, nodes);
            mesh.gather(e, displacement, u);
            double k[solidshell::kDof * solidshell::kDof];
            solidshell::elementStiffness(nodes, steel, solidshell::Formulation::SolidShell, k);
            for (int i = 0; i < solidshell::kDof; ++i) {
                double sum = 0;
                for (int j = 0; j < solidshell::kDof; ++j) sum += k[i * solidshell::kDof + j] * u[j];
                energy += 0.5 * u[i] * sum;
            }
        }
        std::printf("     %2d x %2d elements: %.6e J against %.6e, ratio %.5f\n", refine,
                    8 * refine, energy, want, energy / want);
        last = energy;
    }
    expectNear("the seam's rotational stiffness is 16 D / b per unit length", last, want,
               0.01 * want);
    // The point of the number: it is the *plating's*, and a stiffener on that seam
    // does not change it. Recorded as a limit of the formulation rather than as a
    // success -- a real web has its own torsional and warping stiffness and this
    // has neither.
}

// --- 5. The uniaxial fibre, against its own closed form ---------------------------

void testTheFibreYieldsOnItsFlowCurve() {
    std::printf("\n   one fibre, pulled: the uniaxial return map against the flow curve\n");
    plasticity::Material material = plasticity::shipSteel();
    material.flow = plasticity::linearHardening(355.0e6, 2.0e9);
    material.youngsModulus = 206.0e9;

    // One fibre, one metre long, tied to a pair one metre apart in x and 12 mm in
    // z. Its area is a round number so the force is checkable by inspection.
    const double thickness = 0.012, area = 2.0e-3, span = 1.0;
    constraint::Stiffening stiffening;
    stiffening.material = ah36Steel();
    constraint::Fiber fibre;
    fibre.area = area;
    fibre.offset = 0.0;
    fibre.end[0] = {0, 1, 0.5};
    fibre.end[1] = {2, 3, 0.5};
    stiffening.fiber.push_back(fibre);
    std::vector<double> rest{0, 0, -0.5 * thickness, 0,    0, 0.5 * thickness,
                             span, 0, -0.5 * thickness, span, 0, 0.5 * thickness};
    const constraint::RestFibers forms = constraint::restFibers(stiffening, rest);
    expectNear("the fibre's rest length is the span", forms.length[0], span, 1e-15);

    const auto pull = [&](double strain, int steps, std::vector<constraint::FiberState>& state) {
        double stress = 0;
        for (int s = 1; s <= steps; ++s) {
            std::vector<double> current = rest;
            const double reached = strain * s / steps;
            current[6] += reached * span;
            current[9] += reached * span;
            std::vector<double> force(4 * 3, 0.0);
            constraint::fiberForces(stiffening, forms, current, material, &state, force);
            // The axial force is what the fibre puts on its far end, distributed
            // half and half over that pair by the tie.
            stress = -(force[6] + force[9]) / area;
        }
        return stress;
    };

    // Elastic: sigma = E eps, exactly.
    {
        std::vector<constraint::FiberState> state(1);
        const double strain = 1.0e-4;
        expectNear("below yield the fibre is E eps", pull(strain, 1, state),
                   material.youngsModulus * strain, 1e-9 * material.youngsModulus * strain);
        expectEqual("and nothing yielded", state[0].equivalentPlasticStrain == 0.0 ? 1LL : 0LL,
                    1LL);
    }
    // Plastic: sigma = sigma_y0 + H eps_p with eps_p = eps - sigma/E, which for a
    // linear curve solves to (E sigma_y0 + E H eps) / (E + H).
    const double strain = 0.01;
    const double closed = (material.youngsModulus * 355.0e6 + material.youngsModulus * 2.0e9 * strain) /
                          (material.youngsModulus + 2.0e9);
    std::vector<constraint::FiberState> oneStep(1), manySteps(1);
    const double got1 = pull(strain, 1, oneStep);
    const double got100 = pull(strain, 100, manySteps);
    std::printf("     at 1%% strain: one step %.4f MPa, a hundred steps %.4f MPa, closed form"
                " %.4f MPa\n", got1 / 1e6, got100 / 1e6, closed / 1e6);
    expectNear("the fibre's flow stress is its curve's", got1, closed, 1e-9 * closed);
    // **Step independence**, which is what says the map is solved rather than
    // stepped -- the same property `test_plasticity.cpp` spends its strongest
    // assertion on. A forward-Euler update would drift with the step count.
    expectTrue("and one step reaches it identically to a hundred",
               std::abs(got1 - got100) < 1e-9 * closed);
    expectNear("with the same plastic strain", oneStep[0].equivalentPlasticStrain,
               manySteps[0].equivalentPlasticStrain,
               1e-9 * manySteps[0].equivalentPlasticStrain);
    expectTrue("and it really yielded", oneStep[0].equivalentPlasticStrain > 0.005);

    // Unloading is elastic at E, which is the statement that the plastic strain was
    // stored rather than folded into a secant modulus.
    {
        std::vector<double> current = rest;
        const double back = strain - 1.0e-4;
        current[6] += back * span;
        current[9] += back * span;
        std::vector<double> force(4 * 3, 0.0);
        constraint::fiberForces(stiffening, forms, current, material, &oneStep, force);
        const double unloaded = -(force[6] + force[9]) / area;
        expectNear("unloading falls at E", got1 - unloaded, material.youngsModulus * 1.0e-4,
                   1e-6 * material.youngsModulus * 1.0e-4);
    }

    // **Compression, and then a reversal.** A stiffener on the far side of a dent
    // is in compression, and a fibre that accumulated its plastic strain unsigned
    // would harden identically here and then be wrong on the way back -- the
    // stored plastic strain is what the *next* trial stress is measured from, so
    // the error is invisible until the load turns round. Mutation testing found
    // this exact edit surviving the whole suite.
    {
        std::vector<constraint::FiberState> state(1);
        const double squashed = pull(-strain, 1, state);
        std::printf("     at -1%% strain: %.4f MPa against the closed form's %.4f, with"
                    " eps_p %+.6f\n", squashed / 1e6, -closed / 1e6, state[0].plasticStrain);
        expectNear("compression yields at minus the same flow stress", squashed, -closed,
                   1e-9 * closed);
        expectTrue("and the plastic strain it stored is negative",
                   state[0].plasticStrain < -0.005);
        expectTrue("while the hardening it accumulated is not",
                   state[0].equivalentPlasticStrain > 0.005);
        // The reversal is what an unsigned plastic strain gets wrong: from -1%
        // plastic, coming back to zero *total* strain leaves the fibre in tension
        // at `E * |eps_p|`, capped by the flow stress it has hardened to.
        std::vector<double> current = rest;
        std::vector<double> force(4 * 3, 0.0);
        constraint::fiberForces(stiffening, forms, current, material, &state, force);
        const double reversed = -(force[6] + force[9]) / area;
        const double yielded = plasticity::flowStress(material.flow,
                                                      state[0].equivalentPlasticStrain);
        std::printf("     released to zero total strain it stands at %+.4f MPa, and its flow"
                    " stress is now %.4f MPa\n", reversed / 1e6, yielded / 1e6);
        expectTrue("releasing a compressed fibre leaves it in tension", reversed > 0);
        expectNear("at E times the plastic strain it stored, or at its flow stress if lower",
                   reversed, std::min(-material.youngsModulus * state[0].plasticStrain, yielded),
                   1e-6 * yielded);
    }
}

// A run of exactly two stations is one segment and one fibre per profile station,
// and it has to build: a member that crosses a single element of the patch is a
// member. The check is on `addStiffener` directly because a zone big enough to
// mesh is never that small, and mutation testing found the off-by-one there
// surviving everything.
void testAOneSegmentRunStillBuildsFibres() {
    std::printf("\n   the shortest run a member can have\n");
    const double thickness = 0.012, span = 0.4;
    std::vector<double> rest{0,    0, -0.5 * thickness, 0,    0, 0.5 * thickness,
                             span, 0, -0.5 * thickness, span, 0, 0.5 * thickness};
    constraint::SeamRun run;
    run.bottom = {0, 2};
    run.top = {1, 3};
    run.sign = -1.0;
    constraint::Stiffening stiffening;
    stiffening.material = ah36Steel();
    const std::size_t added =
        constraint::addStiffener(run, flatBar(0.200, 0.010), thickness, rest, stiffening);
    std::printf("     two stations, %.2f m apart: %zu fibres, %.4f m of member, %.4f kg\n", span,
                added, stiffening.length, stiffening.mass);
    expectEqual("two stations build one segment's worth of fibres",
                static_cast<long long>(added), 2LL);
    expectNear("covering the segment", stiffening.length, span, 1e-15);
    expectNear("and weighing the profile's own mass", stiffening.mass,
               profileSection(flatBar(0.200, 0.010)).area * span * ah36Steel().density,
               1e-9 * profileSection(flatBar(0.200, 0.010)).area * span * ah36Steel().density);

    // One station is not a segment and builds nothing, which is the other side of
    // the same boundary.
    constraint::SeamRun single;
    single.bottom = {0};
    single.top = {1};
    constraint::Stiffening nothing;
    nothing.material = ah36Steel();
    expectEqual("one station builds nothing",
                static_cast<long long>(constraint::addStiffener(single, flatBar(0.200, 0.010),
                                                                thickness, rest, nothing)),
                0LL);
}

void testAFibreCarriesNoForceUnderRigidMotion() {
    std::printf("\n   a fibre on a patch that translates and turns\n");
    const StiffenerProfile bar = flatBar(0.200, 0.010);
    const zone::Patch patch = stiffenedPatch(1.4, 0.7, 0.012, 2, 2, bar, true, 2);
    expectTrue("the patch carries fibres", patch.stiffening.fiberCount() > 0);
    const constraint::RestFibers forms =
        constraint::restFibers(patch.stiffening, patch.mesh.position);

    const double angle = 0.4;
    const Vec3 axis = normalize(Vec3{0.2, 0.5, -0.84});
    const Vec3 shift{3.0, -1.0, 2.0};
    std::vector<double> moved(patch.mesh.position.size(), 0.0);
    for (std::size_t node = 0; node < patch.nodeCount(); ++node) {
        const Vec3 p{patch.mesh.position[node * 3], patch.mesh.position[node * 3 + 1],
                     patch.mesh.position[node * 3 + 2]};
        const double c = std::cos(angle), s = std::sin(angle);
        const Vec3 turned = p * c + cross(axis, p) * s + axis * (dot(axis, p) * (1.0 - c)) + shift;
        for (int k = 0; k < 3; ++k) moved[node * 3 + static_cast<std::size_t>(k)] = turned[k];
    }
    std::vector<double> force(patch.nodeCount() * 3, 0.0);
    const constraint::FiberForces rigid = constraint::fiberForces(
        patch.stiffening, forms, moved, plasticity::shipSteel(), nullptr, force);
    double worst = 0;
    for (double f : force) worst = std::max(worst, std::abs(f));
    // Non-vacuous: the same fibres under a stretch do pull, and hard.
    std::vector<double> stretched = patch.mesh.position;
    for (std::size_t node = 0; node < patch.nodeCount(); ++node) stretched[node * 3] *= 1.001;
    std::vector<double> pull(patch.nodeCount() * 3, 0.0);
    const constraint::FiberForces strained = constraint::fiberForces(
        patch.stiffening, forms, stretched, plasticity::shipSteel(), nullptr, pull);
    double pulled = 0;
    for (double f : pull) pulled = std::max(pulled, std::abs(f));

    std::printf("     0.40 rad and 3.7 m: largest nodal force %.3e N, stored energy %.3e J;"
                " a 0.1%% stretch of the same fibres pulls %.3e N\n", worst, rigid.strainEnergy,
                pulled);
    // Not zero, and the reason is arithmetic rather than physics: the fibre's
    // length is a square root of coordinates the rigid motion has moved to 3.7 m,
    // and the tie amplifies their last bit by `1 + 2|w| = 28`, so the strain
    // carries 1e-13 of noise and the fibre's `EA/L` of 8e8 N/m turns that into
    // 1e-4 N. It scales with the *distance translated*, which physics would not.
    // The scale-free statement is the one asserted: it is nine orders below what a
    // tenth of a per cent of real stretch produces.
    expectTrue("a finite rigid motion strains no fibre", worst < 1e-8 * pulled);
    expectTrue("and stores no energy", rigid.strainEnergy < 1e-15);
    expectTrue("but a tenth of a per cent of stretch does", pulled > 1.0e5);
    expectTrue("and stores energy", strained.strainEnergy > 1.0);
}

// --- 6. Mass and the stable step --------------------------------------------------

void testTheFibresMassAndTheStepTheyCost() {
    std::printf("\n   the fibres' mass, and what the tie does to the stable step\n");
    const StiffenerProfile bar = flatBar(0.200, 0.010);
    const StructuralMaterial steel = ah36Steel();
    const double thickness = 0.012, lengthX = 2.0, spanY = 0.7;

    const zone::Patch patch = stiffenedPatch(lengthX, spanY, thickness, 4, 2, bar, true, 2);
    const constraint::RestFibers forms =
        constraint::restFibers(patch.stiffening, patch.mesh.position);
    std::vector<double> mass(patch.nodeCount(), 0.0);
    constraint::lumpFiberMass(patch.stiffening, forms, steel.density, mass);
    double total = 0;
    double negative = 0;
    for (double m : mass) {
        total += m;
        negative = std::min(negative, m);
    }
    const double want = profileSection(bar).area * lengthX * steel.density;
    std::printf("     %.6f kg of stiffener lumped onto the plating, want %.6f, most negative"
                " nodal mass %.3e kg\n", total, want, negative);
    expectNear("the lumping keeps the stiffener's mass exactly", total, want, 1e-12 * want);
    // The reason the lumping is equal-over-the-pair rather than the consistent
    // `T^T M T`: the consistent one has a negative entry for any eccentric tie, and
    // an explicit scheme cannot integrate that. This is the assertion that says the
    // code did not quietly go back to it.
    expectTrue("and every nodal mass is positive", negative == 0.0);

    // The largest eigenvalue of M^-1 K for one fibre, exactly, against a power
    // iteration on the same 12 x 12 block -- an independent route to a number the
    // implementation gets from a closed form.
    std::vector<double> nodal(patch.nodeCount(), 1.7);  // a mass that is not one
    const constraint::FiberStiffness block = constraint::fiberStiffness(
        patch.stiffening.fiber.front(), patch.mesh.position, forms.length.front(),
        steel.youngsModulus);
    const double closed = constraint::fiberFrequencySquared(block, nodal);
    double x[12] = {0.31, -0.17, 0.44, 0.09, -0.63, 0.21, 0.55, 0.02, -0.38, 0.27, 0.14, -0.49};
    double iterated = 0;
    for (int sweep = 0; sweep < 200; ++sweep) {
        double projection = 0;
        for (int i = 0; i < 12; ++i) projection += block.vector[i] * x[i];
        double norm = 0;
        for (int i = 0; i < 12; ++i) {
            x[i] = block.scale * block.vector[i] * projection / nodal[block.dof[i] / 3];
            norm = std::max(norm, std::abs(x[i]));
        }
        iterated = norm;
        for (double& v : x) v /= norm;
    }
    std::printf("     one fibre's largest eigenvalue: closed form %.8e, power iteration %.8e\n",
                closed, iterated);
    expectNear("the rank-one eigenvalue is exact", closed, iterated, 1e-9 * closed);
    expectTrue("and it is a real frequency", closed > 0);

    // What the tie costs the stable step, measured at three element sizes. The
    // plating's own step is thickness governed and flat in the in-plane size
    // (`zone.hpp` §1); a fibre's is not, because EA/L grows as the elements shrink
    // while the mass it is lumped onto falls. So the stiffener re-introduces an
    // in-plane length scale into a step that did not have one, and past some
    // refinement it is the stiffener that sets the step.
    std::printf("     %-12s %14s %14s %8s\n", "subdivision", "plating (s)", "with fibres", "ratio");
    double finestRatio = 1.0;
    for (int subdivision : {2, 4, 8}) {
        const zone::Patch refined =
            stiffenedPatch(lengthX, spanY, thickness, 4, 2, bar, true, subdivision);
        const double ratio = refined.platingTimestep / refined.criticalTimestep;
        std::printf("     %-12d %14.6e %14.6e %8.3f\n", subdivision, refined.platingTimestep,
                    refined.criticalTimestep, ratio);
        expectTrue("the fibres never lengthen the step", refined.criticalTimestep <=
                                                             refined.platingTimestep * (1.0 + 1e-12));
        expectTrue("and the step is usable", refined.criticalTimestep > 0);
        finestRatio = ratio;
    }
    expectTrue("at a fine enough mesh the stiffener sets the step, not the plating",
               finestRatio > 1.5);
    // Non-vacuous the other way: at the coarse mesh the reference ferry actually
    // uses, it does not, or this would be a cost nobody could afford.
    const zone::Patch coarse = stiffenedPatch(lengthX, spanY, thickness, 4, 2, bar, true, 2);
    expectNear("at the reference resolution the plating still sets it",
               coarse.criticalTimestep, coarse.platingTimestep, 1e-12 * coarse.platingTimestep);
}

// --- 7. The mesher's own decisions ------------------------------------------------

void testTheMesherRefusesWhatTheTieCannotCarry() {
    std::printf("\n   what buildPatch does with a member it cannot tie\n");
    const StiffenerProfile bar = flatBar(0.200, 0.010);

    // A web that leans out of the plating's thickness direction has no single
    // eccentricity. Projecting it would be a plausible wrong second moment, so it
    // is refused and said.
    StructuralMesh strip = flatStrip(1.4, 0.7, 0.012, 2, 2, bar, true);
    strip.members[0].rise = normalize(Vec3{0, 1, -1});  // 45 degrees out of the normal
    zone::MeshParams params;
    params.radius = 1000.0;
    params.subdivision = 2;
    params.stiffeners = zone::Stiffeners::Modelled;
    params.outward = {0, 0, 1};
    params.edge = zone::Edge::Free;
    const zone::Patch oblique = zone::buildPatch(strip, {0, 0, 0}, params);
    expectTrue("an oblique web builds no fibres", oblique.stiffening.empty());
    bool said = false;
    for (const std::string& problem : oblique.problems)
        if (problem.find("out of the plating's thickness direction") != std::string::npos)
            said = true;
    expectTrue("and the patch says so", said);

    // A member that is there but out of reach is not silently absent either.
    StructuralMesh far = flatStrip(1.4, 0.7, 0.012, 2, 2, bar, false);
    const zone::Patch none = zone::buildPatch(far, {0, 0, 0}, params);
    bool warned = false;
    for (const std::string& problem : none.problems)
        if (problem.find("no stiffener reached this zone") != std::string::npos) warned = true;
    expectTrue("a zone with no stiffener in it says Modelled did nothing", warned);

    // And the ordinary case still works, so the two checks above are not passing
    // because everything is refused.
    const zone::Patch good = zone::buildPatch(flatStrip(1.4, 0.7, 0.012, 2, 2, bar, true),
                                              {0, 0, 0}, params);
    expectTrue("a perpendicular web does build", !good.stiffening.empty());
    expectEqual("with two fibres per flat bar per segment",
                static_cast<long long>(good.stiffening.fiberCount()),
                2LL * static_cast<long long>(2 * params.subdivision));

    // A member longer than the patch stops where the plating does. The fibres are
    // built over element edges, so a member that leaves the patch and comes back
    // would be two runs rather than one fibre spanning the gap -- and a member
    // that runs off the end simply ends.
    zone::MeshParams narrow = params;
    narrow.radius = 0.5;
    const zone::Patch clipped = zone::buildPatch(flatStrip(2.8, 0.7, 0.012, 4, 2, bar, true),
                                                 {0, 0, 0}, narrow);
    std::printf("     a 0.5 m zone of a 2.8 m strip meshes %.4f m^2 of plating and %.4f m of the"
                " 2.8 m member\n", clipped.area, clipped.stiffening.length);
    expectTrue("the fibres stop where the plating does",
               clipped.stiffening.length > 0 && clipped.stiffening.length < 2.8);
    expectNear("and cover exactly the seam that was meshed", clipped.stiffening.length,
               clipped.area / 0.7, 1e-9 * clipped.area / 0.7);
}

// A member is built over the patch's *element edges*, so one that leaves the
// plating and comes back is two runs and not one fibre spanning the gap. Nothing
// in the ordinary fixtures produces that, and mutation testing found the check
// surviving its removal -- so the geometry is built on purpose: a U of panels
// whose two arms both carry the seam and whose spine does not.
void testAMemberBrokenByThePatchIsTwoRunsAndNotOneFibre() {
    std::printf("\n   a member the patch carries in two pieces\n");
    const double thickness = 0.012, bay = 0.7, half = 0.35;
    StructuralMesh mesh;
    mesh.materials.push_back(ah36Steel());
    mesh.frameSpacing = bay;
    const auto panel = [&](double x0, double x1, double y0, double y1) {
        PlatePanel p;
        p.corner[0] = {x0, y0, 0};
        p.corner[1] = {x1, y0, 0};
        p.corner[2] = {x1, y1, 0};
        p.corner[3] = {x0, y1, 0};
        p.thickness = thickness;
        p.role = PanelRole::Shell;
        mesh.panels.push_back(p);
    };
    // The spine, three bays of it, from y = half to y = 2 * half.
    for (int i = 0; i < 3; ++i) panel(bay * i, bay * (i + 1), half, 2.0 * half);
    // And two arms hanging off its ends, down to y = 0. The middle bay has none,
    // so the seam at y = 0 exists over the arms and not between them.
    panel(0.0, bay, 0.0, half);
    panel(2.0 * bay, 3.0 * bay, 0.0, half);

    StructuralMember member;
    member.a = {0.0, 0.0, 0.0};
    member.b = {3.0 * bay, 0.0, 0.0};
    member.rise = {0, 0, -1};
    member.profile = flatBar(0.200, 0.010);
    member.attachedPlateThickness = thickness;
    member.role = MemberRole::Longitudinal;
    mesh.members.push_back(member);

    zone::MeshParams params;
    params.radius = 1000.0;
    params.subdivision = 2;
    params.stiffeners = zone::Stiffeners::Modelled;
    params.outward = {0, 0, 1};
    params.edge = zone::Edge::Free;
    const zone::Patch patch = zone::buildPatch(mesh, {1.5 * bay, 1.5 * half, 0}, params);
    expectEqualCount("the whole U meshed", patch.panels.size(), std::size_t{5});
    std::printf("     a %.1f m member over a U whose arms cover %.1f m of it: %.4f m of fibres,"
                " %zu of them\n", 3.0 * bay, 2.0 * bay, patch.stiffening.length,
                patch.stiffening.fiberCount());
    // Two arms of one bay each, and *not* the middle bay, which has no plating on
    // the seam at all. A fibre bridging it would be a stiffener hanging in space,
    // and the length is the only place it shows.
    expectNear("the fibres cover the two arms and not the gap between them",
               patch.stiffening.length, 2.0 * bay, 1e-9 * bay);
    expectEqualCount("which is two fibres per element edge over four edges",
                     patch.stiffening.fiberCount(),
                     std::size_t{2 * 2 * static_cast<std::size_t>(params.subdivision)});
    // Two runs, one member. A count that said two would make `members` a property
    // of the zone's shape rather than of the ship's.
    expectEqual("and it is still one member, however many pieces the patch has of it",
                static_cast<long long>(patch.stiffening.members), 1LL);
    // Non-vacuous: the member really does pass through the gap, so a formulation
    // that did not check would have had something to bridge.
    expectTrue("and the member really spans the gap", member.length() > 2.5 * bay);
}

// `solidshell::DofBlock` is the general half of this machinery -- the shape an
// interface coupling between two substructures has as much as a stiffener does --
// and its one failure mode is silent: `BandedSpd::add` drops an entry outside the
// band without a word, so a block coupling distant degrees of freedom would
// simply not be there and the answer would read as slightly soft.
void testADistantDofBlockIsNotDroppedOnTheFloor() {
    std::printf("\n   an extra stiffness block that reaches outside the element bandwidth\n");
    const StructuralMaterial steel = ah36Steel();
    const double length = 2.0, width = 0.4, thickness = 0.02, load = -2.0e4;
    solidshell::HexMesh mesh = solidshell::makePlateMesh(length, width, thickness, 8, 2, 1);
    for (std::size_t node = 0; node < mesh.nodeCount(); ++node)
        if (std::abs(mesh.position[node * 3]) < 1e-9)
            for (int k = 0; k < 3; ++k) mesh.pin(node, k, 0.0);

    // One node at the tip and one at mid-span, both on the +zeta face. They are
    // three x-stations apart, which in this numbering is far outside the element
    // bandwidth -- that is the point.
    std::size_t tip = 0, middle = 0;
    for (std::size_t node = 0; node < mesh.nodeCount(); ++node) {
        const double x = mesh.position[node * 3], y = mesh.position[node * 3 + 1];
        const double z = mesh.position[node * 3 + 2];
        if (std::abs(y) > 1e-9 || z < 0) continue;
        if (std::abs(x - length) < 1e-9) tip = node;
        if (std::abs(x - 0.5 * length) < 1e-9) middle = node;
    }
    expectTrue("the two nodes were found", tip != 0 && middle != 0 && tip != middle);

    std::vector<double> tipLoad(mesh.nodeCount() * 3, 0.0);
    for (std::size_t node = 0; node < mesh.nodeCount(); ++node)
        if (std::abs(mesh.position[node * 3] - length) < 1e-9)
            tipLoad[node * 3 + 2] = load / 6.0;

    std::string problem;
    std::vector<double> free;
    expectTrue("the cantilever solves: " + problem,
               solidshell::solveStatic(mesh, steel, solidshell::Formulation::SolidShell, tipLoad,
                                       free, &problem));
    const double apart = free[tip * 3 + 2] - free[middle * 3 + 2];

    // A stiff spring between the two z degrees of freedom. As its stiffness runs
    // away the two have to agree -- that is a closed-form statement about a
    // penalty tie and needs no second solver to check. An entry dropped for being
    // out of band would leave them exactly where they were.
    solidshell::DofBlock spring;
    const double stiffness = 1.0e12;
    spring.dof = {static_cast<std::uint32_t>(tip * 3 + 2),
                  static_cast<std::uint32_t>(middle * 3 + 2)};
    spring.stiffness = {stiffness, -stiffness, -stiffness, stiffness};
    std::vector<double> tied;
    expectTrue("the same with the block solves: " + problem,
               solidshell::solveStatic(mesh, steel, solidshell::Formulation::SolidShell, {spring},
                                       tipLoad, tied, &problem));
    const double together = tied[tip * 3 + 2] - tied[middle * 3 + 2];
    std::printf("     nodes %zu and %zu, %zu degrees of freedom apart: free %.6e / %.6e m,"
                " tied %.6e / %.6e m (%.3e m apart)\n", middle, tip, 3 * (tip - middle),
                free[middle * 3 + 2], free[tip * 3 + 2], tied[middle * 3 + 2],
                tied[tip * 3 + 2], together);
    expectTrue("the two nodes really move differently without the block",
               std::abs(apart) > 1e-4);
    expectTrue("and the block ties them, so it was assembled and not dropped",
               std::abs(together) < 1e-6 * std::abs(apart));
    // **And it tied them to each other rather than to ground.** `BandedSpd::add`
    // drops an out-of-band entry silently, and a spring whose *off-diagonal* is
    // dropped while its diagonal survives is two enormous springs to earth: both
    // nodes go to zero, they agree perfectly, and the check above passes while the
    // coupling is gone. That is exactly what the mutant that removed the extra
    // blocks' bandwidth contribution did, and it survived until this line existed.
    // The discriminator is *what sets the tied deflection*. A spring between two
    // nodes leaves it set by the plate's own compliance; two springs to earth set
    // it at `load / k`, which here is 2e-8 m and has nothing to do with the plate.
    // Measured: 5.4e-2 m, which is 2.7 million times `load / k`.
    const double toGround = std::abs(load) / stiffness;
    std::printf("     the tied deflection is %.3e m, and a spring to earth would put it at"
                " load/k = %.3e m\n", std::abs(tied[tip * 3 + 2]), toGround);
    expectTrue("the plate's compliance sets the tied deflection, not the spring's stiffness",
               std::abs(tied[tip * 3 + 2]) > 1000.0 * toGround);
    // The coupling has to be outside what the elements alone would have banded, or
    // this passes on a block the element assembly had already made room for.
    expectTrue("and it reached outside the element bandwidth",
               3 * (tip - middle) > 3 * (2 * (2 + 1) + 2));
}

}  // namespace

void runConstraintTests() {
    std::printf("\n--- multi-point constraints and eccentric stiffeners ---\n");
    testTheTieIsExactUnderRigidMotionAndUnderRotationOfThePlate();
    testTheFibresReproduceTheProfileSection();
    testTheStiffenedPanelMatchesItsSectionProperties();
    testTrippingIsNotAZeroEnergyMode();
    testTheSeamRotationStiffnessIsThePlatings();
    testTheFibreYieldsOnItsFlowCurve();
    testAOneSegmentRunStillBuildsFibres();
    testAFibreCarriesNoForceUnderRigidMotion();
    testTheFibresMassAndTheStepTheyCost();
    testTheMesherRefusesWhatTheTieCannotCarry();
    testAMemberBrokenByThePatchIsTwoRunsAndNotOneFibre();
    testADistantDofBlockIsNotDroppedOnTheFloor();
}
