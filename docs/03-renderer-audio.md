# 03 — Renderer and Audio

## Hardware target

Development hardware is a **GTX 1070 Ti (Pascal, 8 GB)**. This is a deliberate
constraint rather than an accident, and it rules things out:

- **No mesh shaders** (Turing+). Geometry pipeline must be compute-driven
  indirect draws, not `VK_EXT_mesh_shader`.
- **No hardware ray tracing.** Lighting is rasterisation plus screen-space and
  probe-based techniques. RT is an optional path added later for hardware that
  has it — not a dependency.
- **No `VK_KHR_dynamic_rendering` performance assumptions** beyond what Pascal
  drivers do well; it is supported but render passes still matter for tiled
  parts of the frame.
- 8 GB VRAM with a ship interior that must be fully resident.

Everything below is chosen to run at 60 fps at 1440p on that card, and to scale up
rather than requiring a scale-down.

## Vulkan 1.3 backend

- **Bindless**: one large descriptor set of textures/buffers, indices passed
  through draw data. Removes per-draw descriptor churn, which is where a
  ship-interior scene with thousands of unique materials would otherwise die.
- **GPU-driven culling**: the full scene lives in GPU buffers; a compute pass does
  frustum, occlusion (HZB two-phase) and small-triangle rejection, and emits
  `vkCmdDrawIndexedIndirectCount`. The CPU submits a handful of draws per frame
  regardless of scene complexity — essential when the CPU is busy solving FEM.
- **Render graph** with automatic barrier and layout transition derivation from
  declared resource usage. Hand-written barriers in a renderer this size are a
  standing bug source.
- **Timeline semaphores** for async compute overlap: the ocean FFT, particle
  simulation, and FEM element pass all run on the compute queue alongside
  graphics.

**Shading**: clustered forward+ with a depth prepass. Chosen over deferred because
of transparency volume (water surfaces, glass, smoke, spray) and because MSAA
stays available, which matters for the thin geometry — railings, ladders, rigging
— that ships are made of.

**Lighting**: physically based, with a strong bias toward the two lighting
situations that dominate: harsh outdoor sun/sky over a highly specular sea, and
enclosed interiors lit by emergency lighting through smoke. Irradiance volumes +
screen-space GI for interiors, analytic sky model outdoors.

## Ocean

- **Tessendorf FFT** spectrum synthesis, several cascades at different spatial
  scales to cover chop through swell without tiling artefacts. Generated on the
  compute queue at 30 Hz.
- **Adaptive projected grid** for the surface, so screen-space triangle density is
  roughly uniform and the horizon does not cost anything.
- **The physics and the visuals read the same spectrum.** Non-negotiable. If the
  wave the ship responds to is not the wave under the bow, the entire premise
  collapses.
- Foam from wave steepness and breaking criterion, advected as a decaying field;
  spray particles from breaking crests and from bow impact; wake from the hull
  (Kelvin pattern plus a simulated near-field).
- Screen-space refraction and reflection, with a planar reflection fallback and
  distance-based cubemap. Subsurface scattering tinting for the shallow water
  colour.
- Underwater: volumetric fog with depth-dependent extinction per wavelength,
  caustics from the surface, and the surface seen from below.

### Ocean — what exists

`engine/gpu/ocean.{hpp,cpp}`, `engine/gpu/shaders/ocean.{vert,frag}`,
`tests/test_ocean.cpp`. A displaced grid patch with a vertex/fragment pipeline
beside `OffscreenRenderer` — no mesh shaders, no ray tracing, nothing Pascal
cannot run.

**It is driven by `sim::WaveField` itself.** Not a visual copy, not a re-seeded
spectrum: `OceanSurface::build` walks the same `WaveComponent` array
`WaveField::elevation()` walks. That is the non-negotiable point above, and it is
asserted rather than asserted-about — at the vertices the grid agrees with
`elevation()` to **1e-6 m**, which is float32 storage and nothing else, and
through the whole render path (project a known world point with `sim::clipToPixel`,
read the elevation the frame actually drew there) to **0.6 mm RMS, 2.2 mm worst**
over 400 samples of a 3 m sea. Normals are the spectrum's own analytic slope, so
the lighting answers to the waves rather than to the tessellation.

**Cost, measured on the shipped path** (one core, `-O2`, 576 components —
Hs 3 m, Tp 9 s, N = 48, M = 12):

| grid | vertices | per frame | per vertex-component |
|---|---|---|---|
| 64 × 64   | 4 225  | 6.5 ms  | 2.67 ns |
| 128 × 128 | 16 641 | 23.6 ms | 2.46 ns |
| 256 × 256 | 66 049 | 92.8 ms | 2.44 ns |

Linear in vertices × components, which is all it is. Two things were learned
getting there:

- **The obvious loop costs 19.8 ns per vertex-component, an eighth of the speed.**
  Along a grid row only x moves, so the phase advances by a constant and
  `(cos ψ, sin ψ)` can be stepped by a fixed rotation instead of re-evaluated —
  one `sincos` per row per component rather than one per vertex, and the sine that
  falls out is exactly what the slope needs, so the normal is free. It is an
  algebraic rearrangement of the same sum: measured against direct evaluation over
  a 257-wide row the worst disagreement is 2e-14 m. Interleaving four independent
  recurrences measured a further 1.6–2.2× and is *not* shipped; it is the next
  thing to do if this path survives.
- **One evaluation per unique grid vertex, not per triangle corner.** Six corners
  reference each interior vertex, so the per-corner version of this costs exactly
  six times as much for the same picture — the same redundancy the physics tick
  was caught paying. `tests/test_ocean.cpp` asserts the vertex and index counts so
  a regression to per-corner shows up as a failure rather than as a slow frame.

**Resolution is a correctness constraint, not a quality knob**, and the criterion
is not the one you would reach for. `engine/core/geometry.hpp` records the same
hazard on the hull side. Measured at Hs 3 m, Tp 9 s, N = 48 — dominant wavelength
126 m, **shortest component 11.4 m** — worst cell-centre error against the
analytic field:

| cells per dominant wavelength | 16 | 32 | 64 | 128 |
|---|---|---|---|---|
| grid spacing | 7.90 m | 3.95 m | 1.98 m | 0.99 m |
| worst error, % of Hs | 24% | 8.2% | 2.2% | 0.56% |

**Resolving the dominant wavelength is not enough** — sixteen cells across it
still invents a quarter of the wave height, because the short components are below
that grid's Nyquist. The criterion is the shortest component the discretisation
carries, and about eight cells across it is where the surface stops inventing
height: 1.4 m spacing for that sea, 80 000 vertices over a 400 m patch, and
**92 ms per frame on one core**. `shortestWavelength()` exists so callers compute
that rather than guess it.

Which is the argument for the FFT cascade above, now with a number on it. Note
that evaluating the spectrum in a *vertex shader* is not the fix: it is the same
sum with the same components, and it re-does the per-vertex `sincos` the row
recurrence just removed. The FFT is the fix, and it still reads this spectrum.

Deliberately absent: **specular and Fresnel.** The shading is Lambert against a
directional sun plus a hemispheric sky term and is view-independent, which is what
makes a flat sea exactly one colour over its whole surface — a closed form the
tests check to the last bit, and the strongest single assertion in the file
because it drives the entire path against an answer statable in advance. Specular
is what makes water look like water and it costs that assertion, so it arrives
with the reflection work rather than ahead of it.

### Ocean — what the tests are pointed at

`tests/test_ocean.cpp`. A grid of displaced quads looks like the sea whatever it
is doing, so nothing is eyeballed: a flat sea's exact colour and the shoelace area
of its projected quadrilateral, the sign flip of a wave train after half a period,
the amplitude a long-crested sea's crest must reach, and the linear-interpolation
error `a (1 − cos(k h / 2)) cos ψ`, which is an identity rather than a bound and is
asserted both for a grid that resolves the wave and for one that does not.

Three of them earned their place:

- **The load-bearing check's tolerance was derived from the artefact under test.**
  It bounded the discrepancy partly by measuring the grid's own interpolation error
  with `OceanSurface::sampleElevation` — so a surface with a systematic elevation
  error inflated its own tolerance by exactly the amount it was wrong by. A
  mutation settled it: scaling every elevation by 1.02 was caught by four other
  assertions in the file and sailed straight through this one. The bound now comes
  only from the wave field and the camera, the sub-pixel offset of the pixel centre
  is *solved for* by Newton rather than swallowed by tolerance, and the same
  mutation is now rejected by 374 of 400 samples. Worst bound fell from 31 mm to
  6 mm on the way.
- **A comparison that cannot fail proves nothing**, so the same check runs against
  the field one second later and must reject nearly every sample. Mutation-tested
  the same way: a renderer 50 ms behind the physics trips ten assertions here.
- **The first vertex-agreement check was measuring the wrong quantity.** It
  compared against `elevation()` evaluated at the vertex's stored float32 position,
  while the builder had evaluated the field at the unrounded grid point — so it was
  measuring position quantisation times local slope, and failed at 1.6e-6 m against
  a 1e-6 m expectation. The expectation was right and the comparison was wrong;
  position is now asserted separately, against its own float32 bound.

## Hull deformation and damage

The FEM produces displaced node positions; the renderer must show them without a
CPU round-trip.

- Visual mesh is skinned to the FEM mesh with precomputed weights; deformation is
  applied in a compute shader directly from the solver's node buffer, which never
  leaves the GPU.
- Tears change topology. Handled by a compute pass that rebuilds the affected
  index ranges into a preallocated slack region — no reallocation mid-frame.
- Torn edges get procedural jagged geometry and exposed-metal material blending
  driven by plastic strain, so a dent looks like a dent and a rupture looks like
  torn steel.
- Paint, rust and scorch are separate texture-space layers updated by the sim
  (corrosion state, temperature history), not baked variants.

## Water in compartments

The interior water is the hardest rendering problem in the project, because it is
seen from *inside* at close range in enclosed spaces.

- Quiescent compartments: an analytic plane clipped to the room, which is exactly
  what the sim already computes. Cheap and perfectly stable.
- Dynamic compartments: screen-space fluid surface reconstruction from particles
  (depth buffer → bilateral smoothing → normal reconstruction), plus thickness for
  absorption. The classic Green/van der Laan approach; it holds up at close range
  and costs a few passes.
- The transition between the two must be invisible. This constrains the demotion
  criteria in `02-simulation.md` §6 as much as the physics does.
- Wetness: surfaces above the current and historical waterline get a wetness mask
  that darkens albedo and sharpens specular, decaying over time.

## Fire and smoke

Volumetric raymarching against the gas solver's density and temperature fields,
with blackbody emission from the temperature directly — so the colour of a fire is
a physical consequence of how hot it is, not an art choice. Sparse volume
representation (NanoVDB-style) so only the burning compartments cost anything.
Smoke layer descent in the multi-zone model renders as a genuine stratified layer,
which is both correct and much cheaper than a full volume.

## VR

Supported, and it constrains the above:

- Stereo via multiview (`VK_KHR_multiview`), single-pass.
- 90 Hz floor with a hard budget; the effects above have explicit VR variants
  (screen-space fluid at reduced resolution, volumetrics at half-res with
  temporal reprojection).
- Motion sickness is a genuine design problem in a game about a rolling deck.
  Mitigations: comfort options for horizon-locked reference frames, vignetting
  under acceleration, and — because the ship's motion is *simulated* rather than
  animated — the option to lock the camera to the world horizon while the deck
  moves under you, which is what a real person's vestibular system does anyway.
- OpenXR, not a vendor SDK.

## Audio

Ships are a rich acoustic environment and it is mostly structure-borne.

- **Sources from the sim, not from triggers**: hull plating groans from FEM stress
  levels; the flooding roar's intensity from actual volumetric flow rate through
  each opening; engine note from actual RPM and load; the specific sound of water
  finding a new path when a door fails.
- **Propagation**: a graph over the compartment network (the same graph the water
  and gas use), with per-edge transmission loss depending on whether the door is
  open, closed or under water. Sound getting quieter as a bulkhead closes is free
  if the topology is already there.
- Reverb from compartment volume and material absorption, computed rather than
  authored per room.
- HRTF binaural rendering, and underwater filtering with the correct speed-of-sound
  change.
- Structure-borne transmission: an impact forward is heard aft through the steel
  before it is heard through the air.
