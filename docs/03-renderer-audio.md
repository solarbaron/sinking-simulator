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

## Hull rendering and materials

### What exists

`engine/gpu/hull.{hpp,cpp}`, `engine/gpu/material.{hpp,cpp}`,
`engine/gpu/shaders/hull.{vert,frag}`, `engine/gpu/materials/marine.materials`,
`tests/test_hull_render.cpp`. A lit solid: `sim::Ship`'s hull shaded by an
analytic metallic-roughness BRDF, on a plain vertex/fragment pipeline over indexed
triangles with one storage buffer for the material table. Nothing Pascal cannot
run.

**The hull and the sea go through one pipeline**, and that is the design rather
than an optimisation. `SceneMesh` is the geometry container for everything the lit
pass draws; `appendOcean` converts an `OceanSurface` into it with a material index
of its own. The ship is therefore *in* the water rather than composited beside it,
the occlusion between them is the depth test doing its job, and the sea is a
material like any other — which makes "add a material without recompiling" true of
the water too. This supersedes the flat-shaded path `tools/ferry_view` drives
through `OffscreenRenderer`; `MeshVertex` carries a baked colour and no normal, so
there is nothing there for a light to answer to.

Triangles are drawn in index order, so appending the sea *after* the hull is what
the tests do deliberately: a broken depth test then paints the water straight over
the ship instead of hiding behind the draw order.

### The BRDF

Metallic-roughness, which is what glTF, Filament and every DCC tool already speak,
so an imported material means something here. Per channel, for normal `n`, unit
vector `l` toward the light, unit vector `v` toward the eye and `h = normalize(l + v)`:

```
L = base * sky * (0.5 + 0.5 n.z)  +  (base (1 - metal) / pi  +  D * Vis * F) * max(n.l, 0) * sun

alpha = roughness^2
D     = a2 / (pi ((n.h)^2 (a2 - 1) + 1)^2)                    GGX / Trowbridge-Reitz
Vis   = 0.5 / (n.l sqrt((n.v)^2 (1 - a2) + a2) + n.v sqrt((n.l)^2 (1 - a2) + a2))
                                                              Smith, height-correlated
F     = F0 + (1 - F0) (1 - v.h)^5,   F0 = mix(0.04, base, metal)     Schlick
```

and the stored pixel is `round(255 clamp(exposure * L, 0, 1))`.

Three choices, each of which buys a closed form and each of which is a trade rather
than an oversight:

- **The diffuse lobe is exactly Lambert** — no `(1 - F)` factor, which glTF's
  reference BRDF carries and Filament drops. Keeping it couples the diffuse term
  to the view direction and costs the exact cosine law the tests sweep over
  nineteen angles. It is worth a few per cent at grazing incidence and the sweep is
  worth more.
- **No tone mapping and no sRGB curve.** The pixel is a clamped linear radiance, so
  every assertion is an equality rather than an inversion of a curve. Same trade the
  sea's specular records: the transfer function arrives with the post chain.
- **The hemispheric sky term is not scaled by `(1 - metal)`.** A metal has no
  diffuse lobe but does mirror the sky, and its mirror reflectance *is* its base
  colour, so `base * sky` stands in for that until image-based lighting exists.

Roughness below 0.03 is **rejected, not clamped**. At `alpha = 9e-4` the GGX peak is
already 4e5 and the lobe is a fraction of a pixel wide, so a smaller number is not a
smoother surface, it is a sampling problem — and silently clamping a mod's value is
how data stops meaning what it says.

### Materials are data

There is no material table in any header and no enumeration in any shader. The
shipped set is `engine/gpu/materials/marine.materials`:

```
material painted_steel_topside
    base_colour 0.68 0.70 0.72
    roughness   0.35
    metalness   0.00
```

`MaterialLibrary::load` merges any number of such files and the fragment shader
indexes an std430 storage buffer built from whatever was loaded, so **a mod adds a
surface by loading a second file and naming it**. A repeated name overrides in
place, which means an index resolved before the mod loaded still points at the same
surface. `docs/05-data-modding-validation.md` commits to a ship being a directory;
a paint scheme is part of a ship, and this is that commitment kept on the render
side.

Shipped: painted steel topside, boot-topping, antifouling, bare steel (the only
conductor — paint over steel is a dielectric however metallic the plate under it
is), rusted steel (an oxide, so *not* a conductor), timber deck, glass, and the sea.

**The loader fails closed.** An unknown key, a number out of range, a duplicate name
inside one file, or a block that stops before it is complete leaves the library
exactly as it was, with a line-numbered reason. That is the defect CLAUDE.md records
in `World::load`, and it is checked here against every line truncation *and* every
byte truncation of a valid file.

Which material a hull triangle gets is a `HullPaint`: material names and the z bands
they occupy, decided per triangle from its centroid and geometric normal. **Bands are
decided in the body frame, not the world frame** — paint is on the hull, so a ship
heeled thirty degrees still has her boot-topping where it was painted.

### Normals

Flat (the geometric normal, on all three corners) or smooth (area-weighted average
over the faces meeting at a welded position, except across an edge sharper than the
crease angle). Area weighted rather than plain averaged, so a fan of slivers cannot
outvote the one large face it sits against. A ship needs both: a smooth turn of the
bilge and a hard deck edge.

### Cost, measured on a GTX 1070 Ti

The whole test ship (796 hull triangles) plus a 400 m sea patch, one frame,
best of several:

| target | vertices | triangles | CPU build | upload | **GPU** | submit | readback |
|---|---|---|---|---|---|---|---|
| 512 × 384 | 6 613 | 8 988 | 0.07 ms | 0.10 ms | **0.02 ms** | 0.12 ms | 4.2 ms |
| 512 × 384 | 19 029 | 33 564 | 0.26 ms | 0.25 ms | **0.03 ms** | 0.14 ms | 4.5 ms |
| 512 × 384 | 39 637 | 74 524 | 0.54 ms | 0.46 ms | **0.04 ms** | 0.21 ms | 4.2 ms |
| 1920 × 1080 | 39 637 | 74 524 | — | 0.47 ms | **0.26 ms** | 0.95 ms | 12.1 ms |

**The GPU is not where the frame goes, and at this scale nothing is.** 0.26 ms at
1080p is 1.6% of a 60 Hz frame for a whole ship and a whole sea, so none of this
needs optimising before Phase 3 lands on top of it. Two things the table is
actually saying:

- **Readback dominates the offscreen path and is not a frame cost.** It is a
  `vkCmdCopyImageToBuffer` plus a host memcpy of the whole colour attachment, which
  exists so tests can assert on pixels. An on-screen renderer does none of it.
- **The remaining CPU cost is the sea, not the hull.** `OceanSurface::build` (which
  the ocean section already measures) plus converting its vertices into the shared
  format; the hull's own build is under 0.05 ms and does not move. The conversion is
  a straight copy and is the price of one pipeline for both, which is worth paying
  for a shared depth buffer and is the first thing to remove if the sea ever moves
  to a compute-generated buffer.

Cull mode is `NONE`, for the reason the sea's is: a hull is seen from inside once it
is cut away or flooded, and a back face there is a real surface. The fragment shader
flips the normal for the far side. On a closed hull the front faces win the depth
test anyway, so this costs fill rate and no correctness — and at these numbers the
fill rate is not worth the class of handedness bug that culling invites.

### Hull rendering — what the tests are pointed at

`tests/test_hull_render.cpp`. A lit hull looks like a ship whatever it is doing —
more so than the sea does, because the silhouette carries the recognition on its own
and the shading can be arbitrarily wrong underneath it. So the BRDF above is written
out a second time in the test file, from the formula rather than from the shader,
and every shaded assertion is a pixel predicted before the render:

- a hull face pointing **exactly** at the sun, and one **exactly** ninety degrees off
  it, which is `base * sky * 0.5` and nothing else — one multiplication, and the
  strongest single line in the file;
- **Lambert's cosine law over nineteen angles**, with the specular lobe *cancelled
  analytically* rather than modelled: the lobe is achromatic for a dielectric, so the
  red-minus-blue difference is pure diffuse and
  `(red - blue) / (base_r - base_b) - sky` is exactly `(sun / pi) cos(theta)`.
  Measured worst discrepancy 0.006 over the sweep; a `cos^2` falloff is rejected by
  15 of the 20 samples and a linear one by 17, because "darker as the light swings
  away" is satisfied by both;
- the material a known hull point is painted with, found through `sim::clipToPixel`
  and confirmed against the depth channel so an occluded sample is skipped rather
  than compared against a different triangle;
- the far surface drawn **second** and still losing the depth test;
- port and starboard views of a symmetric hull mirroring — with the sun constrained
  to the centreplane, because otherwise the two sides are genuinely lit differently
  and the images have no business agreeing.

Two readback channels make the geometric claims integers rather than colour
comparisons: **`MaterialId`** (which material won the depth test, 16-bit) and
**`Depth`** (clip-space z, 16-bit). The occlusion check is then exact — a sea point
beyond the ship must come back as hull *and* at a depth nearer than the water's,
and a sea point in front must come back as water at the depth the camera matrix
predicts.

**The instrument that is not a pixel comparison is an energy balance.** A BRDF that
reflects more light than arrives is unphysical, and no pixel comparison can see it
when the shader and the prediction make the same mistake — which is the pattern
CLAUDE.md's whole table is about. Integrating the specular lobe over the hemisphere
with `F0` driven to one gives a directional albedo bounded above by one:

| roughness | 0.20 | 0.35 | 0.50 | 0.70 | 1.00 |
|---|---|---|---|---|---|
| albedo at n·v = 0.60 | 0.997 | 0.965 | 0.872 | 0.695 | 0.412 |

The shortfall is the single-scattering energy loss GGX is known for and is the
multi-scatter term's job; the first version of this check demanded 0.45 as a floor
and failed at 0.32, which was the **expectation** being wrong rather than the model.
Helmholtz reciprocity is asserted alongside it, to 1e-12.

### What mutation testing changed

Twelve deliberate defects were introduced and the suite re-run. Most were caught
loudly — a diffuse lobe missing its `1/pi` trips 14 assertions, a flipped depth
compare 46, a material index off by one 33, paint decided in the world frame exactly
1. **Five escaped entirely**, and each one is now a test that did not exist before:

| escaped | now caught by |
|---|---|
| smooth normals not area weighted | a ridge with one side eight times the area of the other and the same triangle count on each — every earlier construction had equal-area triangles at a shared vertex, so weighting made no difference |
| back faces not flipped | a plate viewed from behind; nothing else here looks at a hull from inside |
| Schlick's exponent 5 → 4 | a grazing view at n·v = 0.13, where it is worth 11 codes; every other probe sat near the normal, where the two agree to a rounding error |
| the sea's indices not rebased onto its own vertices | a structural check on the appended index range. Un-rebased they still address *valid* vertices — the hull's — so they draw a plausible surface, and against a **flat calm sea they draw the identical surface**, because shuffling the corners of a plane gives back the plane. The calm sea that makes the occlusion predictions exact is precisely what blinds it to this |
| the sea's analytic normal replaced by flat up | a bit-identical comparison against `OceanSurface`'s own normals, with a guard that they actually tilt |

That last one is the general lesson repeated: the configuration that makes an
assertion exact is often the configuration that makes it blind, and the two have to
be separated deliberately.

### Deliberately absent

- **Transparency.** `Material::opacity` is carried because it is a property of the
  material and the data must not lie about the ship, but the pass ignores it: there
  is no sort and no refraction, so glass renders as the smooth dark dielectric it
  is from outside. This is the same shape of decision as the sea's missing specular.
- **Shadows, and any indirect light beyond the hemispheric term.** A ship casts a
  hard shadow on her own deck and it is not here.
- **Sky reflection.** Which is why `sea_water`'s base colour is deliberately above
  the physical two or three per cent: at the real value the sea renders very nearly
  black, correctly and uselessly. It comes back down when the reflection work lands,
  and being a number in a file rather than a constant in a shader is the point of it
  living there.
- **A straight paint-band boundary between two of the hull's own waterlines.**
  Assignment is per triangle, so a band edge that does not fall on a waterline of the
  offsets table is only as straight as the tessellation — a sawtooth, at exactly the
  line a person looks at hardest. Fairing band edges onto the waterlines is what a
  real hull is painted to anyway and it is what the test ship does; the general fix
  is the texture-space paint layer this document already plans for rust and scorch.

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
