// SPDX-License-Identifier: MIT
#include "promotion.hpp"

#include "solid_shell.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>

namespace sim::promotion {
namespace {

using solidshell::kDof;
using solidshell::kNodes;

// Robust rather than mean, because the flat-ship guard is the whole point: one
// enormous station must not drag the background up behind it. Empty is zero,
// which makes every excess test pass -- correct, since with no stations there is
// no ship to be flat.
double median(std::vector<double> values) {
    if (values.empty()) return 0.0;
    const std::size_t half = values.size() / 2;
    std::nth_element(values.begin(), values.begin() + static_cast<std::ptrdiff_t>(half),
                     values.end());
    return values[half];
}

// Linear interpolation of the bending moment between stations; the ends clamp.
double momentAt(const HullGirder& girder, double x) {
    if (girder.stations.empty()) return 0.0;
    if (x <= girder.stations.front().x) return girder.stations.front().moment;
    if (x >= girder.stations.back().x) return girder.stations.back().moment;
    for (std::size_t i = 1; i < girder.stations.size(); ++i) {
        const GirderStation& lo = girder.stations[i - 1];
        const GirderStation& hi = girder.stations[i];
        if (x > hi.x) continue;
        const double span = hi.x - lo.x;
        if (!(span > 0)) return lo.moment;
        const double t = (x - lo.x) / span;
        return lo.moment + t * (hi.moment - lo.moment);
    }
    return girder.stations.back().moment;
}

// The panel whose centroid is nearest a point. -1 when there are none.
int nearestPanel(const StructuralMesh& structure, const Vec3& point) {
    int best = -1;
    double nearest = 0;
    for (std::size_t i = 0; i < structure.panels.size(); ++i) {
        const double d = length(structure.panels[i].centroid() - point);
        if (best < 0 || d < nearest) {
            nearest = d;
            best = static_cast<int>(i);
        }
    }
    return best;
}

// Utilisation of each trigger at the station nearest `x`, for the hold test --
// which is deliberately geometry-free, so keeping a zone costs nothing to decide.
struct StationUtilisation {
    double yield = 0, buckle = 0, collapse = 0;
};

template <class T>
const T* nearestStation(const std::vector<T>& stations, double x) {
    const T* best = nullptr;
    double nearest = 0;
    for (const T& s : stations) {
        const double d = std::abs(s.x - x);
        if (best == nullptr || d < nearest) {
            nearest = d;
            best = &s;
        }
    }
    return best;
}

StationUtilisation utilisationAt(const TierZero& tier, double x) {
    StationUtilisation out;
    if (const GirderStress* s = nearestStation(tier.stress, x)) out.yield = s->utilisation;
    if (const GirderBuckling* b = nearestStation(tier.buckling, x)) out.buckle = b->utilisation;
    if (const StrengthStation* c = nearestStation(tier.strength, x))
        out.collapse = c->margin > 0 ? 1.0 / c->margin : 0.0;
    return out;
}

double costOf(int elements, double thickness, const Criterion& criterion) {
    const double scale = thickness > 0 ? criterion.costThickness / thickness : 1.0;
    return criterion.coreSecondsPerElement * static_cast<double>(elements) * scale;
}

}  // namespace

const char* name(Trigger trigger) {
    switch (trigger) {
        case Trigger::None: return "none";
        case Trigger::Yield: return "yield";
        case Trigger::Buckling: return "buckling";
        case Trigger::Collapse: return "collapse";
        case Trigger::Contact: return "contact";
    }
    return "?";
}

// --- What Tier-0 knows ---------------------------------------------------------

TierZero tierZero(const Ship& ship, const Sea& sea, const StructuralMesh& structure,
                  const Scantlings& scantlings, const TierZeroParams& params) {
    const auto begin = std::chrono::steady_clock::now();
    TierZero out;

    out.yieldStrength = params.yieldStrength > 0 ? params.yieldStrength
                        : !structure.materials.empty() ? structure.materials.front().yieldStrength
                                                       : ah36Steel().yieldStrength;
    if (!(out.yieldStrength > 0)) {
        out.problems.push_back("no yield strength: every utilisation would be a ratio to zero");
        return out;
    }

    out.girder = hullGirder(ship, sea, params.stations);
    for (std::string& problem : validateGirder(out.girder)) out.problems.push_back(problem);

    out.stress = girderStress(out.girder, structure, out.yieldStrength);
    out.buckling = girderBuckling(out.stress, structure, scantlings);
    if (params.collapse)
        out.strength = longitudinalStrength(out.girder, structure, scantlings, params.shedExponent,
                                            params.curvatureSteps);

    out.yieldUtilisation = worstUtilisation(out.stress, &out.yieldX);
    out.buckleUtilisation = worstBucklingUtilisation(out.buckling, &out.buckleX);
    for (const StrengthStation& s : out.strength) {
        if (!(s.margin > 0)) continue;
        const double utilisation = 1.0 / s.margin;
        if (utilisation > out.collapseUtilisation) {
            out.collapseUtilisation = utilisation;
            out.collapseX = s.x;
        }
    }

    const auto end = std::chrono::steady_clock::now();
    out.seconds = std::chrono::duration<double>(end - begin).count();
    return out;
}

// --- Pieces of the criterion ---------------------------------------------------

double panelSpan(const PlatePanel& panel) {
    // Each opposite pair averaged, so a tapered panel reports its mean side rather
    // than its narrow end -- and the *shorter* pair, because that is the direction
    // the plating between stiffeners actually spans.
    const double a = 0.5 * (length(panel.corner[1] - panel.corner[0]) +
                            length(panel.corner[2] - panel.corner[3]));
    const double b = 0.5 * (length(panel.corner[3] - panel.corner[0]) +
                            length(panel.corner[2] - panel.corner[1]));
    return std::min(a, b);
}

double platingCollapsePressure(double yieldStrength, double thickness, double span) {
    if (!(span > 0) || !(thickness > 0) || !(yieldStrength > 0)) return 0.0;
    const double ratio = thickness / span;
    return 4.0 * yieldStrength * ratio * ratio;
}

int estimateElements(const StructuralMesh& structure, const Vec3& impact,
                     const zone::MeshParams& mesh) {
    if (!(mesh.radius > 0) || mesh.subdivision < 1) return 0;
    int panels = 0;
    for (const PlatePanel& panel : structure.panels) {
        if (panel.role != mesh.role) continue;
        if (length(panel.centroid() - impact) > mesh.radius) continue;
        ++panels;
    }
    return panels * mesh.subdivision * mesh.subdivision;
}

std::vector<Candidate> candidates(const StructuralMesh& structure, const TierZero& tier,
                                  const std::vector<ContactPatch>& contacts,
                                  const Criterion& criterion) {
    std::vector<Candidate> out;

    // The three girder triggers share a shape: a per-station utilisation, a median
    // over the stations that had one, an absolute threshold and a local excess.
    // Only the fibre they point at differs, so they are driven from one loop.
    struct Girder {
        double x = 0, utilisation = 0;
        bool deckFibre = false;
        Trigger trigger = Trigger::None;
    };
    std::vector<Girder> sites;
    std::vector<double> yieldValues, buckleValues, collapseValues;

    for (const GirderStress& s : tier.stress) {
        yieldValues.push_back(s.utilisation);
        Girder g;
        g.x = s.x;
        g.utilisation = s.utilisation;
        g.deckFibre = std::abs(s.stressDeck) >= std::abs(s.stressKeel);
        g.trigger = Trigger::Yield;
        sites.push_back(g);
    }
    for (const GirderBuckling& b : tier.buckling) {
        buckleValues.push_back(b.utilisation);
        Girder g;
        g.x = b.x;
        g.utilisation = b.utilisation;
        g.deckFibre = b.deckInCompression;
        g.trigger = Trigger::Buckling;
        sites.push_back(g);
    }
    for (const StrengthStation& s : tier.strength) {
        if (!(s.margin > 0)) continue;
        collapseValues.push_back(1.0 / s.margin);
        Girder g;
        g.x = s.x;
        g.utilisation = 1.0 / s.margin;
        // Collapse starts on the compressed fibre, and the moment says which:
        // hogging arches the hull and compresses the keel.
        g.deckFibre = s.appliedMoment < 0;
        g.trigger = Trigger::Collapse;
        sites.push_back(g);
    }

    const double yieldBackground = median(yieldValues);
    const double buckleBackground = median(buckleValues);
    const double collapseBackground = median(collapseValues);

    for (const Girder& g : sites) {
        double promote = 0, background = 0;
        switch (g.trigger) {
            case Trigger::Yield:
                promote = criterion.yieldPromote;
                background = yieldBackground;
                break;
            case Trigger::Buckling:
                promote = criterion.bucklePromote;
                background = buckleBackground;
                break;
            case Trigger::Collapse:
                promote = criterion.collapsePromote;
                background = collapseBackground;
                break;
            default: continue;
        }
        if (!(promote > 0)) continue;
        if (g.utilisation < promote) continue;
        if (g.utilisation - background < criterion.localExcess) continue;

        // Only now is the section worth 0.09 ms: this is the one place the girth
        // position of a zone gets decided, and it is decided from a beam, which is
        // the criterion's second stated miss.
        const HullGirderSection section = hullGirderSection(structure, g.x);
        if (!(section.modulusDeck > 0) || !(section.modulusKeel > 0)) continue;
        const Vec3 target{g.x, 0.0, g.deckFibre ? section.zDeck : section.zKeel};
        const int panel = nearestPanel(structure, target);
        if (panel < 0) continue;

        Candidate c;
        c.panel = panel;
        c.impact = structure.panels[static_cast<std::size_t>(panel)].centroid();
        c.x = g.x;
        c.role = structure.panels[static_cast<std::size_t>(panel)].role;
        c.trigger = g.trigger;
        c.utilisation = g.utilisation;
        c.background = background;
        c.score = g.utilisation / promote;
        c.why = std::string(name(g.trigger)) + " " + std::to_string(g.utilisation) + " against a " +
                std::to_string(background) + " background on the " +
                (g.deckFibre ? "deck" : "keel") + " fibre";
        out.push_back(c);
    }

    for (const ContactPatch& contact : contacts) {
        if (!(contact.force > 0) || !(contact.radius > 0)) continue;
        const int panel = nearestPanel(structure, contact.centre);
        if (panel < 0) continue;
        const PlatePanel& struck = structure.panels[static_cast<std::size_t>(panel)];
        const double yieldStrength =
            static_cast<std::size_t>(struck.material) < structure.materials.size()
                ? structure.materials[static_cast<std::size_t>(struck.material)].yieldStrength
                : tier.yieldStrength;
        const double collapse =
            platingCollapsePressure(yieldStrength, struck.thickness, panelSpan(struck));
        if (!(collapse > 0)) continue;
        const double pressure = contact.force / (kPi * contact.radius * contact.radius);
        const double utilisation = pressure / collapse;
        if (!(criterion.contactPressure > 0)) continue;
        if (utilisation < criterion.contactPressure) continue;

        Candidate c;
        c.panel = panel;
        c.impact = contact.centre;
        c.x = contact.centre.x;
        c.role = struck.role;
        c.trigger = Trigger::Contact;
        c.utilisation = utilisation;
        c.background = 0;  // a beam has no background for a load it cannot represent
        c.score = utilisation / criterion.contactPressure;
        c.why = "contact at " + std::to_string(pressure / 1e6) + " MPa against a bay that hinges at " +
                std::to_string(collapse / 1e6) + " MPa";
        out.push_back(c);
    }

    // One zone per panel: keep the strongest claim on it. Deterministic ties, by
    // panel index, so two runs of the same load promote the same patches.
    std::sort(out.begin(), out.end(), [](const Candidate& a, const Candidate& b) {
        if (a.panel != b.panel) return a.panel < b.panel;
        if (a.score != b.score) return a.score > b.score;
        return static_cast<int>(a.trigger) < static_cast<int>(b.trigger);
    });
    out.erase(std::unique(out.begin(), out.end(),
                          [](const Candidate& a, const Candidate& b) { return a.panel == b.panel; }),
              out.end());

    for (Candidate& c : out) {
        zone::MeshParams mesh = criterion.mesh;
        mesh.role = c.role;
        c.elements = estimateElements(structure, c.impact, mesh);
        c.cost = costOf(c.elements, structure.panels[static_cast<std::size_t>(c.panel)].thickness,
                        criterion);
    }

    std::sort(out.begin(), out.end(), [](const Candidate& a, const Candidate& b) {
        if (a.score != b.score) return a.score > b.score;
        return a.panel < b.panel;
    });
    return out;
}

// --- The state machine ---------------------------------------------------------

Promoter::Promoter(Criterion criterion) : criterion_(criterion) {}

zone::MeshParams Promoter::meshFor(const Candidate& candidate) const {
    zone::MeshParams mesh = criterion_.mesh;
    mesh.role = candidate.role;
    return mesh;
}

int Promoter::activeElements() const {
    int total = 0;
    for (const Active& a : active_) total += a.elements;
    return total;
}

double Promoter::activeCost() const {
    double total = 0;
    for (const Active& a : active_) total += a.cost;
    return total;
}

void Promoter::clear() {
    active_.clear();
    qualifying_.clear();
}

Review Promoter::review(const StructuralMesh& structure, const TierZero& tier,
                        const std::vector<ContactPatch>& contacts) {
    const auto begin = std::chrono::steady_clock::now();
    Review out;
    ++reviews_;
    out.considered = candidates(structure, tier, contacts, criterion_);

    const auto isCandidate = [&](int panel) {
        for (const Candidate& c : out.considered)
            if (c.panel == panel) return true;
        return false;
    };

    // Does an *already promoted* zone still deserve to be? The hold thresholds,
    // and no local-excess test: that is the hysteresis, and it is what stops a
    // zone that has just raised the background from demoting itself.
    const auto holds = [&](const Active& zone) {
        if (isCandidate(zone.panel)) return true;
        const StationUtilisation u = utilisationAt(tier, zone.x);
        if (u.yield >= criterion_.yieldHold && criterion_.yieldHold > 0) return true;
        if (u.buckle >= criterion_.buckleHold && criterion_.buckleHold > 0) return true;
        if (u.collapse >= criterion_.collapseHold && criterion_.collapseHold > 0) return true;
        for (const ContactPatch& contact : contacts) {
            if (!(contact.force > 0) || !(contact.radius > 0)) continue;
            if (length(contact.centre - zone.impact) > criterion_.mesh.radius) continue;
            const PlatePanel& struck = structure.panels[static_cast<std::size_t>(zone.panel)];
            const double yieldStrength =
                static_cast<std::size_t>(struck.material) < structure.materials.size()
                    ? structure.materials[static_cast<std::size_t>(struck.material)].yieldStrength
                    : tier.yieldStrength;
            const double collapse =
                platingCollapsePressure(yieldStrength, struck.thickness, panelSpan(struck));
            if (!(collapse > 0)) continue;
            const double pressure = contact.force / (kPi * contact.radius * contact.radius);
            if (pressure / collapse >= criterion_.contactHold) return true;
        }
        return false;
    };

    std::vector<Active> kept;
    for (Active& zone : active_) {
        if (holds(zone)) {
            zone.idleReviews = 0;
            kept.push_back(zone);
            continue;
        }
        ++zone.idleReviews;
        if (zone.idleReviews >= std::max(1, criterion_.hold)) {
            out.demoted.push_back(zone);
            ++demotions_;
        } else {
            kept.push_back(zone);
        }
    }
    active_ = kept;

    // Dwell: consecutive reviews a candidate has qualified for. A candidate that
    // drops out loses its count entirely rather than decaying, because the whole
    // point is *consecutive*.
    std::vector<std::pair<int, int>> nextQualifying;
    for (const Candidate& c : out.considered) {
        int previous = 0;
        for (const auto& entry : qualifying_)
            if (entry.first == c.panel) previous = entry.second;
        nextQualifying.emplace_back(c.panel, previous + 1);
    }
    std::sort(nextQualifying.begin(), nextQualifying.end());
    qualifying_ = nextQualifying;

    const auto dwellOf = [&](int panel) {
        for (const auto& entry : qualifying_)
            if (entry.first == panel) return entry.second;
        return 0;
    };
    const auto alreadyActive = [&](int panel) {
        for (const Active& a : active_)
            if (a.panel == panel) return true;
        return false;
    };
    const auto tooClose = [&](const Vec3& point) {
        const double limit = criterion_.separation * criterion_.mesh.radius;
        for (const Active& a : active_)
            if (length(a.impact - point) < limit) return true;
        return false;
    };

    int budgetUsed = activeElements();
    bool refused = false;
    for (const Candidate& c : out.considered) {
        if (alreadyActive(c.panel)) continue;
        if (dwellOf(c.panel) < std::max(1, criterion_.dwell)) continue;
        if (tooClose(c.impact)) continue;
        if (criterion_.elementBudget > 0 && budgetUsed + c.elements > criterion_.elementBudget) {
            refused = true;
            continue;
        }
        Active zone;
        zone.panel = c.panel;
        zone.impact = c.impact;
        zone.x = c.x;
        zone.role = c.role;
        zone.trigger = c.trigger;
        zone.score = c.score;
        zone.elements = c.elements;
        zone.cost = c.cost;
        zone.promotedAtReview = reviews_;
        active_.push_back(zone);
        out.promoted.push_back(zone);
        ++promotions_;
        budgetUsed += c.elements;
    }
    if (refused)
        out.problems.push_back("the element budget refused a qualifying zone; nothing already"
                               " running was evicted for it");

    out.elementsActive = activeElements();
    out.costActive = activeCost();
    const auto end = std::chrono::steady_clock::now();
    out.microseconds = std::chrono::duration<double>(end - begin).count() * 1e6;
    return out;
}

// --- The pre-load --------------------------------------------------------------

PreloadCheck preloadFor(const HullGirder& girder, const StructuralMesh& structure,
                        const zone::Patch& patch, double obliquityLimit) {
    PreloadCheck out;
    if (patch.empty()) {
        out.problems.push_back("no patch to pre-load");
        return out;
    }
    const HullGirderSection section = hullGirderSection(structure, patch.centre.x);
    if (!(section.secondMoment > 0)) {
        out.problems.push_back("no hull girder section at the patch's station");
        return out;
    }
    out.moment = momentAt(girder, patch.centre.x);
    out.neutralAxis = section.neutralAxis;

    zone::Preload preload;
    preload.reference = section.neutralAxis;
    preload.gradient = out.moment / section.secondMoment;
    preload.stress = 0.0;
    out.surfaceStress = preload.at(patch.centre.z);

    // How far the patch's own normal leans out of the athwartships plane. The
    // uniaxial state is traction-free only when it lies in it; what is left on the
    // face otherwise is `sigma sin^2(phi)`, which is a real unbalanced load and
    // would ring rather than sit.
    const double axisLength = length(patch.axis);
    const double along = axisLength > 0 ? std::abs(patch.axis.x) / axisLength : 1.0;
    out.obliquity = std::asin(std::min(1.0, std::max(0.0, along)));
    out.tractionError = std::abs(out.surfaceStress) * along * along;

    if (out.obliquity > obliquityLimit) {
        out.problems.push_back("the patch's normal leans " + std::to_string(out.obliquity) +
                               " rad out of the athwartships plane, so a longitudinal stress is not"
                               " traction-free on it: not pre-loaded");
        return out;
    }
    out.preload = preload;
    out.applied = preload.active();
    return out;
}

// --- The reaction back ---------------------------------------------------------

SectionReduction reactionOf(const StructuralMesh& structure, const zone::Patch& patch,
                            const zone::Solver& solver) {
    SectionReduction out;
    if (patch.empty()) {
        out.problems.push_back("no patch: nothing to reduce");
        return out;
    }

    const std::vector<solidshell::ElementPlasticState>& state = solver.elementState();
    const std::vector<double>& current = solver.position();
    const std::vector<double>& rest = solver.rest();
    const bool haveState = state.size() == patch.elementCount();
    if (!haveState)
        out.problems.push_back("the solve carried no plastic state, so nothing can have torn:"
                               " the reduction is thinning only");

    // Per panel: the carrying thickness its meshed elements have left, area
    // weighted, plus the dent and the torn area for reporting.
    struct Accumulator {
        double carrying = 0;   // sum over elements of area * (t_now / t_rest)
        double meshed = 0;     // sum over elements of area
        double torn = 0;       // m^2
        double thinning = 0;   // m, the worst
        double outOfPlane = 0; // m, the worst
    };
    std::vector<Accumulator> byPanel(patch.panels.size());
    const auto slot = [&](int panel) -> Accumulator* {
        for (std::size_t i = 0; i < patch.panels.size(); ++i)
            if (patch.panels[i] == panel) return &byPanel[i];
        return nullptr;
    };

    for (std::size_t e = 0; e < patch.elementCount(); ++e) {
        Accumulator* acc = slot(patch.panelOf[e]);
        if (acc == nullptr) continue;
        const double area = patch.elementArea[e];
        acc->meshed += area;

        if (haveState && state[e].torn) {
            acc->torn += area;
            continue;  // carries nothing
        }

        double nodes[kDof];
        patch.mesh.gather(e, current, nodes);
        double inPlane = 0, thickness = 0;
        solidshell::elementSize(nodes, &inPlane, &thickness);

        // Against the thickness this element had **when the zone was promoted**,
        // not against the nominal plate thickness. `patch.mesh.position` is exactly
        // that configuration, since a `Preload` moves the rest state and leaves the
        // current one where the mesher put it. They differ on a curved or a
        // distorted element, where volume over mid-surface area is not the
        // thickness the strake was authored with, and that difference would be
        // reported as damage the punch did not do.
        double promotedNodes[kDof];
        patch.mesh.gather(e, patch.mesh.position, promotedNodes);
        double promotedInPlane = 0, promotedThickness = 0;
        solidshell::elementSize(promotedNodes, &promotedInPlane, &promotedThickness);

        // Clamped at one: a section *reduction* that could add material back is
        // not a reduction, and bending does locally thicken an element.
        const double ratio = promotedThickness > 0
                                 ? std::min(1.0, std::max(0.0, thickness / promotedThickness))
                                 : 1.0;
        acc->carrying += area * ratio;
        acc->thinning = std::max(acc->thinning, promotedThickness * (1.0 - ratio));

        double restNodes[kDof];
        patch.mesh.gather(e, rest, restNodes);
        for (int a = 0; a < kNodes; ++a) {
            const Vec3 moved{nodes[a * 3] - restNodes[a * 3], nodes[a * 3 + 1] - restNodes[a * 3 + 1],
                             nodes[a * 3 + 2] - restNodes[a * 3 + 2]};
            acc->outOfPlane = std::max(acc->outOfPlane, std::abs(dot(moved, patch.axis)));
        }
    }

    out.xLo = 1e300;
    out.xHi = -1e300;
    for (std::size_t i = 0; i < patch.panels.size(); ++i) {
        const Accumulator& acc = byPanel[i];
        const int index = patch.panels[i];
        if (index < 0 || static_cast<std::size_t>(index) >= structure.panels.size()) {
            out.problems.push_back("the patch names a panel the structure does not have");
            continue;
        }
        const PlatePanel& panel = structure.panels[static_cast<std::size_t>(index)];
        const double full = panel.area();
        if (!(full > 0)) continue;

        PanelDamage damage;
        damage.panel = index;
        damage.meshedFraction = std::min(1.0, acc.meshed / full);
        damage.outOfPlane = acc.outOfPlane;
        damage.thinning = acc.thinning;
        damage.tornArea = acc.torn;
        // The part of the panel the zone never meshed counts as intact, because
        // nothing looked at it. Anything else would let the zone radius decide how
        // much of the ship is damaged -- the defect `indentation.hpp` records.
        //
        // `min(meshed, full)` and the snap below are not tidiness. A panel meshed
        // edge to edge has `meshed` equal to `full` to rounding, so a *completely*
        // torn one came out at an effectiveness of 1e-16 rather than zero -- and
        // `reduce()` then left a plate 1e-18 m thick, which is not zero, so
        // `collapseElementsAt` kept it, and a plate that thin has a critical stress
        // of 1e-16 Pa. That drove the whole section's first-yield curvature to
        // nothing and reported the ferry's ultimate moment as 1e-21 N m. An
        // arithmetic epsilon read as a ship with no strength at all.
        const double intact = std::max(0.0, full - std::min(acc.meshed, full));
        damage.effectiveness = std::min(1.0, std::max(0.0, (acc.carrying + intact) / full));
        if (damage.effectiveness < 1e-6) damage.effectiveness = 0.0;
        // A part in ten thousand of a 12 mm plate is a micron, which is below the
        // mill tolerance the strake was rolled to. Reporting it as lost section is
        // noise, and it is noise a pre-loaded zone makes *more* of: measured on the
        // ferry, an unfloored reduction named twice as many damaged panels under
        // her own 13 MPa of hogging stress as it did without it, every one of the
        // extra ones intact to six figures.
        if (damage.effectiveness >= 1.0 - 1e-4) continue;

        out.panels.push_back(damage);
        out.worstEffectiveness = std::min(out.worstEffectiveness, damage.effectiveness);
        out.worstOutOfPlane = std::max(out.worstOutOfPlane, damage.outOfPlane);
        out.lostPlateArea += (1.0 - damage.effectiveness) * full;
        const double density =
            static_cast<std::size_t>(panel.material) < structure.materials.size()
                ? structure.materials[static_cast<std::size_t>(panel.material)].density
                : ah36Steel().density;
        out.lostSteelMass += (1.0 - damage.effectiveness) * full * panel.thickness * density;
        for (const Vec3& corner : panel.corner) {
            out.xLo = std::min(out.xLo, corner.x);
            out.xHi = std::max(out.xHi, corner.x);
        }
    }
    if (out.panels.empty()) {
        out.xLo = out.xHi = 0;
        out.worstEffectiveness = 1.0;
    }
    std::sort(out.panels.begin(), out.panels.end(),
              [](const PanelDamage& a, const PanelDamage& b) { return a.panel < b.panel; });
    return out;
}

StructuralMesh reduce(const StructuralMesh& structure, const SectionReduction& reduction) {
    StructuralMesh out = structure;
    for (const PanelDamage& damage : reduction.panels) {
        if (damage.panel < 0 || static_cast<std::size_t>(damage.panel) >= out.panels.size())
            continue;
        const double factor = std::min(1.0, std::max(0.0, damage.effectiveness));
        out.panels[static_cast<std::size_t>(damage.panel)].thickness *= factor;
    }
    return out;
}

double dentedCompressiveCapacity(double yieldStrength, double criticalStress, double deviation,
                                 double thickness) {
    if (!(yieldStrength > 0) || !(criticalStress > 0)) return 0.0;
    if (!(thickness > 0) || !(deviation > 0)) return std::min(yieldStrength, criticalStress);
    // sigma^2 - sigma (sigma_cr (1 + eta) + sigma_y) + sigma_y sigma_cr = 0, the
    // smaller root. At eta = 0 the discriminant is exactly (sigma_cr - sigma_y)^2
    // and the root is min(sigma_cr, sigma_y), which is the continuity that makes
    // this a knockdown rather than a separate model.
    const double eta = 6.0 * deviation / thickness;
    const double b = criticalStress * (1.0 + eta) + yieldStrength;
    const double discriminant = b * b - 4.0 * yieldStrength * criticalStress;
    if (!(discriminant > 0)) return 0.0;
    return 0.5 * (b - std::sqrt(discriminant));
}

}  // namespace sim::promotion
