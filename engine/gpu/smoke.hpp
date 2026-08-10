// SPDX-License-Identifier: MIT
//
// The participating medium a two-zone fire actually has, and nothing else.
//
// `docs/06-roadmap.md` Phase 4's "volumetric fire and smoke rendering", and
// `docs/03-renderer-audio.md`'s "smoke layer descent in the multi-zone model
// renders as a genuine stratified layer, which is both correct and much cheaper
// than a full volume".
//
// --- What the simulation contains, and what it therefore cannot be asked to draw
//
// `engine/sim/fire.hpp` is a **two-zone model**. Per compartment it carries a hot
// upper layer over a cool lower one: two masses, two internal energies, two
// species loadings, and a horizontal interface between them. That is eight
// numbers and a box. It is not a field, and the difference is not one of
// resolution:
//
//   * **There is no plume.** Heskestad's correlation appears in `fire.hpp` as an
//     entrainment *rate* -- one scalar, kg/s -- and as a mean flame height. The
//     plume's shape is not stored, not solved for, and not recoverable: a
//     correlation is a relation between integrals, and the geometry that produced
//     it was integrated away when it was fitted. A rising column drawn here would
//     be an invention, and this file does not draw one.
//   * **There is no ceiling jet, and no horizontal structure of any kind.** A zone
//     is *well mixed by definition*. Smoke does not bank against the deckhead here
//     and does not thin toward a doorway, because the model says it does not.
//   * **There is no vertical structure inside a layer.** Each layer is exactly
//     uniform, so the medium is two homogeneous slabs and the transfer integral is
//     closed form rather than a ray march. That is a consequence of the physics,
//     not an optimisation.
//   * **The interface is a plane, and it is drawn sharp.** `fire.hpp` says the
//     interface is horizontal in the *body* frame and that a heeled compartment's
//     layer geometry is unsolved. So it tilts with the ship and it has no
//     thickness. Softening it would be drawing a mixing layer the model does not
//     have -- and it would cost the assertion in `tests/test_smoke_render.cpp`
//     that the interface lands on the pixel `sim::clipToPixel` predicts.
//   * **There is no flame.** `DesignFire` carries a heat release, a base height
//     and a pan diameter; `Plume::flameHeight` returns Heskestad's *mean* height.
//     A flame is intermittent at that height half the time and the model carries
//     no flame temperature at all. A glowing cylinder of the mean height would be
//     a scalar wearing a shape. What burns, here, is the layer: a grey gas at the
//     temperature the model solved for.
//
// So what is drawn is exactly this: **two homogeneous, emitting, absorbing slabs
// per compartment**, the upper one hot and sooty, the lower one cool and mostly
// clear, separated by a plane at the height the model reports.
//
// --- The prism, and why it is not the compartment's mesh -----------------------
//
// `Model::attach` derives each gas space from the compartment's bounding box in z
// and its air volume, and sets `floorArea = gasVolume / height`. The model's own
// space is therefore a **prism**: uniform plan area between floor and deckhead,
// which is what makes `interfaceZ()` a single number at all. This file draws that
// prism -- the bounding box's aspect ratio in plan, scaled so the plan area is
// exactly `floorArea` -- and so the drawn gas volume equals `gasVolume` and the
// drawn upper-layer volume equals `upperVolume()` to machine precision.
//
// Drawing the compartment's real mesh instead would put gas where the zone model
// never solved for any. Measured on the ferry, the difference is small where it
// matters and large where it does not: the engine rooms fill 99.0% of their
// bounding box and the holds 96.7-98.7%, while the forepeak fills 57.8% and the
// forward wing tanks 58.0%. A machinery-space fire is drawn on essentially the
// right box; a forepeak fire is drawn on the box the *model* used, which is
// smaller than the space, and that is the model's approximation showing rather
// than the renderer's.
//
// --- The optics ----------------------------------------------------------------
//
// One slab of extinction coefficient `k` and path length `d`, at uniform
// temperature, over a background radiance `L_bg`, is the exact solution of the
// radiative transfer equation for a non-scattering grey medium:
//
//     L = B (1 - exp(-k d)) + exp(-k d) L_bg
//
// -- Beer-Lambert transmittance and Kirchhoff emissivity, which are the same
// `exp(-k d)`. **A transparent gas cannot glow**: at `k = 0` the emissivity is
// exactly zero however hot the layer is, so a hot but soot-free layer contributes
// nothing, and the composite is bit-for-bit the background.
//
// Two slabs compose in the order the ray meets them, and along a straight ray
// through a horizontally stratified medium the *order is unambiguous* even when
// the geometry chops the path into pieces: z is monotone in the ray parameter, so
// all of the upper-layer path precedes all of the lower-layer path (or follows
// it), and the pieces are separated only by vacuum, which neither emits nor
// absorbs. That is why two path lengths and a sign are enough, and why there is
// no march.
//
// **Scattering is not modelled and its absence is not free.** Smoke scatters as
// much as it absorbs; treating extinction as pure absorption makes a lit smoke
// layer darker than it should be, and makes it not glow when a torch is shone
// into it. Named rather than hidden -- an in-scattering term needs a light list,
// which is a renderer this one is not yet.
//
// The colour of the emission is Planck's law, integrated over three contiguous
// 100 nm bands -- R [600, 700], G [500, 600], B [400, 500] nm. It is a *spectral
// binning, not a colourimetric transform*: a CIE integration would put the
// chromaticity somewhere slightly different, and there is no CIE data in this
// repository to do it with. Equal-width contiguous bands make the triple a fair
// sample of the visible spectrum, Wien's displacement law holds on it, and the
// whole-spectrum integral recovers Stefan-Boltzmann -- all three are asserted in
// `tests/test_smoke.cpp`.
//
// SI units and the body frame per CLAUDE.md; Vulkan clip space (y down, z in
// [0, 1]) for everything that reaches a pixel.
#pragma once

#include "../core/math.hpp"
#include "../core/png.hpp"
#include "../sim/fire.hpp"
#include "hull.hpp"

#include <cstdint>
#include <string>
#include <vector>

struct VkRenderPass_T;
struct VkFramebuffer_T;
struct VkPipeline_T;
struct VkPipelineLayout_T;
struct VkDescriptorSetLayout_T;
struct VkDescriptorPool_T;
struct VkDescriptorSet_T;
struct VkQueryPool_T;
struct VkSampler_T;

namespace gpu {

// ---------------------------------------------------------------------------
// Blackbody
// ---------------------------------------------------------------------------

// The 2019 SI defining constants, exact by definition. Stefan-Boltzmann is
// *derived* from them rather than quoted, on the same terms as `fire.hpp`'s
// caloric constants: a quoted 5.670374419e-8 against these h, c and k would
// disagree in the ninth digit, and `tests/test_smoke.cpp` asserts the quadrature
// against this sigma at 1e-9, which is tight enough to see that.
inline constexpr double kPlanck = 6.62607015e-34;      // J s
inline constexpr double kLightSpeed = 299792458.0;     // m/s
inline constexpr double kBoltzmann = 1.380649e-23;     // J/K
// sigma = 2 pi^5 k^4 / (15 c^2 h^3).
double stefanBoltzmann();

// Planck's spectral radiance, W / (m^2 sr m). `wavelength` in metres.
double planckRadiance(double wavelength, double temperature);

// Blackbody radiance over the whole spectrum, `sigma T^4 / pi`, W / (m^2 sr).
// The closed form the band quadrature is measured against.
double blackbodyRadiance(double temperature);

// `planckRadiance` integrated over [lo, hi] metres of wavelength, W / (m^2 sr).
//
// Composite 8-point Gauss-Legendre in **log wavelength**, because the spectrum
// spans six decades and a uniform rule over the whole of it under-resolves the
// peak. `intervals` is exposed so a test can run a refinement study rather than
// trust a number.
double blackbodyBandRadiance(double lo, double hi, double temperature, int intervals = 64);

// The three display bands, W / (m^2 sr) each. Contiguous and equal width, so the
// triple is a spectral sample rather than a colourimetric one.
inline constexpr double kBandEdges[4] = {400e-9, 500e-9, 600e-9, 700e-9};
// Intervals `blackbodyBands` uses over each 100 nm band. Four, and the number is
// measured rather than chosen: see `smoke.cpp`, and `tests/test_smoke.cpp` holds
// it against a 4096-interval reference. This is the only expression in the file
// that runs per layer per frame, which is why it is worth being tight about.
inline constexpr int kDisplayBandIntervals = 4;
void blackbodyBands(double temperature, double out[3]);

// ---------------------------------------------------------------------------
// The medium
// ---------------------------------------------------------------------------

// How the species tracer and the layer temperature become optics. One policy for
// the whole model: these are properties of smoke, not of a compartment.
struct SmokeShading {
    // Mass extinction coefficient of the tracer, m^2/kg. The SFPE Handbook's
    // recommended value for the smoke of flaming combustion (Mulholland and
    // Croarkin) is 8700 m^2/kg, and that is a coefficient **for soot**.
    //
    // `fire.hpp`'s tracer is "combustion products", generated at a default yield
    // of 0.05 kg per kg of fuel -- which is a soot-like yield for a sooty fuel and
    // a total-products yield for a clean one. This is the one knob that
    // distinguishes those two readings of the same number, and it is here rather
    // than per compartment because it is a property of what is burning.
    double massExtinction = 8700.0;

    // The gas temperature whose brightest band reads exactly 1.0. **This is a
    // display mapping and the only one in this file**; every ratio below it --
    // between channels, between layers, between temperatures -- is Planck's.
    //
    // It has to exist because the lit pass's radiance scale is arbitrary
    // (`SceneView::sunColour` is about pi for a sun that is 2e7 W/(m^2 sr)), and
    // it has to be a temperature rather than a gain because the exponential is
    // brutal: the visible-band radiance of a grey body rises by a factor of ~500
    // between 900 K and 1200 K, so no linear mapping shows both. 1100 K sits in
    // the middle of a post-flashover upper layer.
    double referenceTemperature = 1100.0;   // K
};

// One well-mixed layer as the renderer sees it: an extinction coefficient and an
// emitted radiance, both already resolved out of the gas state.
struct SmokeLayer {
    double extinction = 0;               // 1/m
    double emission[3] = {0, 0, 0};      // display units, the lit pass's scale
};

// One compartment's gas: the prism the zone model solves on, in the ship's body
// frame, plus the transform that puts it in the world.
struct SmokeVolume {
    std::string name;
    sim::Vec3 lo, hi;            // body frame, the model's own prism
    double interfaceZ = 0;       // body frame; horizontal in *body*, per fire.hpp
    SmokeLayer upper, lower;
    sim::Mat3 rotation;          // body -> world; identity by default (Mat3{} is I)
    sim::Vec3 translation;       // body -> world

    // The plan area of `lo`..`hi`. Equal to the gas model's `floorArea` when this
    // came from `volumesFromFire`, which is what makes the drawn gas volume the
    // model's gas volume.
    double planArea() const { return (hi.x - lo.x) * (hi.y - lo.y); }
    sim::Vec3 centreWorld() const;
};

// The optics of one layer, from the gas state. Split out from the loop below so
// a test can state the chain -- concentration, extinction, emissivity, Planck --
// on one layer and one temperature rather than on a whole ship.
SmokeLayer layerFrom(const sim::fire::Layer& layer, double volume, const SmokeShading& shading);

// The emitted radiance of a grey body at `temperature`, in display units:
// `blackbodyBands(T)` divided by the brightest band of
// `blackbodyBands(referenceTemperature)`. A layer at the reference temperature
// therefore reads exactly 1.0 in its brightest channel.
void emissiveColour(double temperature, const SmokeShading& shading, double out[3]);

// Every tracked gas space as a drawable prism, in the ship's current attitude.
//
// A gas compartment whose `shipCompartment` is `sim::kSea` -- a bare model in a
// test, with no ship behind it -- gets a **square** footprint of the model's own
// floor area centred on the body origin, because a square is the only footprint
// that does not claim to know something. Everything else takes the compartment's
// bounding box aspect ratio in plan.
std::vector<SmokeVolume> volumesFromFire(const sim::fire::Model& model, const sim::Ship& ship,
                                         const SmokeShading& shading);

// ---------------------------------------------------------------------------
// The transfer integral, in closed form
// ---------------------------------------------------------------------------

// Beer-Lambert. `exp(-k d)`, and exactly 1.0 when either argument is zero.
double transmittance(double extinction, double distance);

// Where a ray meets the prism, in the ray's own parameter (metres, because
// `direction` is required to be a unit vector and the body transform is a
// rotation). False when it misses, or when the intersection is entirely behind
// the eye or beyond `maxDistance`.
bool intersectVolume(const SmokeVolume& volume, const sim::Vec3& origin,
                     const sim::Vec3& direction, double maxDistance, double& tEnter,
                     double& tExit);

// How much of [tEnter, tExit] lies in each layer, and which the ray meets first.
//
// `upperFirst` is decided by the sign of the ray's body-frame z component, not by
// where the eye is: a ray that descends meets the upper layer first wherever it
// started. A ray with no z component lies wholly in one layer and the flag is
// then not read.
struct LayerPath {
    double upper = 0;         // m
    double lower = 0;         // m
    bool upperFirst = true;
};
LayerPath layerPath(const SmokeVolume& volume, const sim::Vec3& origin,
                    const sim::Vec3& direction, double tEnter, double tExit);

// The exact two-slab solution, over `background`. The independent statement of
// what `shaders/smoke.frag` computes: `tests/test_smoke_render.cpp` predicts
// every pixel from here and holds the GPU against it, which is the same
// arrangement `material.hpp`'s BRDF has with `hull.frag`.
void compositeOver(const SmokeVolume& volume, const LayerPath& path, const double background[3],
                   double out[3]);


// ---------------------------------------------------------------------------
// The camera, recovered from the matrix rather than passed beside it
// ---------------------------------------------------------------------------

// The view direction and the two constants that turn a depth-buffer value back
// into a distance along it: `viewDistance = b / (depth + a)`.
//
// **Derived from the mvp, not taken from the caller.** A perspective
// `projection * view` has row 3 equal to `(-forward, dot(forward, eye))` and row
// 2 equal to `a` times row 3 plus `(0, 0, 0, b)`, so both constants and the
// direction are already in the matrix. Passing them alongside would be a second
// copy of the camera that could disagree with the first, which is exactly the
// failure `eyeAgreesWithCamera` exists to catch on the eye.
//
// False for an orthographic projection, which has no centre of projection and no
// such row.
struct DepthBasis {
    sim::Vec3 forward;      // unit, world
    double a = 0, b = 0;
};
bool depthBasisFrom(const float mvp[16], DepthBasis& out);

// ---------------------------------------------------------------------------
// The renderer
// ---------------------------------------------------------------------------

// Where a frame went, split at the pass boundary, because the question the
// volumetric pass raises is how much it costs *on top of* the lit solid.
struct SmokeFrameCost {
    double uploadSeconds = 0;
    double opaqueGpuSeconds = 0;   // timestamps around the lit pass
    double smokeGpuSeconds = 0;    // timestamps around the volumetric pass
    double submitSeconds = 0;
    double readbackSeconds = 0;
    double totalSeconds = 0;
    std::size_t vertices = 0;
    std::size_t triangles = 0;
    std::size_t volumes = 0;
};

// The lit solid with a participating medium composited over it.
//
// **Two passes, one command buffer.** The first is `HullRenderer`'s pipeline --
// the same `hull.vert.spv` and `hull.frag.spv`, not a copy of them -- so an empty
// volume list produces a frame that is *bit-identical* to `HullRenderer::render`,
// and `tests/test_smoke_render.cpp` asserts exactly that. The second reads the
// first's depth buffer as a texture to find where each ray stops, and composites
// with premultiplied alpha, which makes `exp(-k d)` the destination blend factor
// and therefore makes Beer-Lambert the hardware's arithmetic rather than the
// shader's.
//
// Volumes are sorted back to front by the distance from the eye to their centres
// and drawn in that order. Emission makes the composite order-dependent even for
// media that do not overlap -- the nearer one's glow is attenuated by less -- so
// this is a correctness requirement rather than a transparency convention. It is
// exact for media that do not overlap along a ray, which the zone model's prisms
// do not, being derived from compartments that a pairwise overlap check already
// keeps apart.
//
// Pascal target: vertex/fragment over a generated box, one storage buffer, one
// sampled depth image. No mesh shaders, no ray tracing, no ray march.
class SmokeRenderer {
public:
    ~SmokeRenderer();

    bool create(Device& device, std::uint32_t width, std::uint32_t height,
                const std::string& shaderDirectory, std::string& error);
    void destroy();
    bool valid() const { return smokePipeline_ != nullptr; }

    std::uint32_t width() const { return width_; }
    std::uint32_t height() const { return height_; }

    // An empty `volumes` is legal and is the bit-identity case above. False when
    // `mvp` is not a perspective projection, because the depth reconstruction has
    // no meaning then and a silently wrong distance is worse than a refusal.
    bool render(const float mvp[16], const SceneView& view, const SceneMesh& mesh,
                const MaterialLibrary& library, const std::vector<SmokeVolume>& volumes,
                const float clearColour[4], core::Image& out);

    const SmokeFrameCost& lastFrame() const { return lastFrame_; }

private:
    bool ensureGeometryCapacity(std::size_t vertexBytes, std::size_t indexBytes);
    bool ensureMaterials(const MaterialLibrary& library);
    bool ensureVolumes(const std::vector<SmokeVolume>& sorted);

    Device*             device_ = nullptr;
    std::uint32_t       width_ = 0;
    std::uint32_t       height_ = 0;

    Device::Image2D     colour_;
    Device::Image2D     depth_;
    VkRenderPass_T*     opaquePass_ = nullptr;
    VkRenderPass_T*     smokePass_ = nullptr;
    VkFramebuffer_T*    opaqueFramebuffer_ = nullptr;
    VkFramebuffer_T*    smokeFramebuffer_ = nullptr;

    VkDescriptorSetLayout_T* materialLayout_ = nullptr;
    VkDescriptorSetLayout_T* volumeLayout_ = nullptr;
    VkDescriptorPool_T* descriptorPool_ = nullptr;
    VkDescriptorSet_T*  materialSet_ = nullptr;
    VkDescriptorSet_T*  volumeSet_ = nullptr;
    VkSampler_T*        depthSampler_ = nullptr;

    VkPipelineLayout_T* opaqueLayout_ = nullptr;
    VkPipeline_T*       opaquePipeline_ = nullptr;
    VkPipelineLayout_T* smokeLayout_ = nullptr;
    VkPipeline_T*       smokePipeline_ = nullptr;
    VkQueryPool_T*      queryPool_ = nullptr;

    Device::Buffer      vertexBuffer_;
    Device::Buffer      indexBuffer_;
    Device::Buffer      materialBuffer_;
    Device::Buffer      volumeBuffer_;
    Device::Buffer      readback_;

    std::uint64_t       materialRevision_ = 0;
    std::size_t         materialCapacity_ = 0;
    std::size_t         volumeCapacity_ = 0;
    SmokeFrameCost      lastFrame_;
    std::vector<SmokeVolume> sorted_;
};

}  // namespace gpu
