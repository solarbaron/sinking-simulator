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
  roughly uniform and the horizon does not cost anything. The geometric cascade
  below is the first half of this: it makes the horizon nearly free, but the
  triangle density is uniform in *world* space rather than in screen space.
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
cannot run — and, over it, a **cascade of concentric rings** that carries the same
patch out to the horizon for a cost that grows with the *logarithm* of the reach.

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

Which is why a *uniform* patch cannot reach the horizon, and why the sea's edge
was visible against the sky in `seaway_view`: applying that criterion out to the
horizon costs the square of the distance, so the patch stopped at 525 m and the
frame showed it. Note that evaluating the spectrum in a *vertex shader* is not
the fix either: it is the same sum with the same components, and it re-does the
per-vertex `sincos` the row recurrence just removed.

### The cascade

`OceanCascade` / `OceanCascadeSurface` in `engine/gpu/ocean.{hpp,cpp}`.
Concentric square rings: level 0 is a full grid of `resolution` cells, and each
level out **doubles both the cell size and the extent**. A ring is therefore the
same cell count minus the quarter its predecessor already covers — three quarters
of level 0's cost for four times the area — and

> **reach is exponential in the number of levels while cost is linear in it.**

That is the whole of the idea. Ten levels reach 512 times as far as one. The
criterion above then applies where it means something and stops applying where it
does not.

**How far is far enough is a fact about the camera, not a number.** A flat sea's
horizon *is* the eye's own horizontal plane, so a surface point at horizontal
distance `D` sits `atan(h / D)` below it, and a patch that stops at `D` leaves
exactly that much sky between its edge and the water that should still be there.
One pixel at the centre of the frame subtends `2 tan(fov / 2) / pixelHeight`, so

```
oceanHorizonReach(h, fov, pixels) = h * pixels / (2 tan(fov / 2))
```

is where the edge disappears into the horizon line. Doubling the eye height
doubles the reach required and so does doubling the vertical resolution; nothing
about the sea state enters it, because a crest is below the eye whatever it is
doing. `tests/test_ocean.cpp` asserts the *consequence* through a real camera
matrix rather than the formula against itself: at that reach the edge projects
**exactly one pixel** below the vanishing point of the sea plane, and at a quarter
of it exactly four. Both are identities, so both are asserted as equalities.

**Cracking is prevented by construction, not by tolerance.** Two levels meeting
at a boundary is the whole difficulty of the technique, and there are two ways for
it to go wrong:

- A **T-junction**: the fine level has a vertex halfway along each of the coarse
  level's boundary edges, and a coarse triangle that interpolates straight past it
  leaves a wedge of background showing through.
- A **step**: two levels evaluating the same point and disagreeing, because they
  carry different components or simply because two floating-point sums of the same
  numbers in a different order are not the same number.

Both are removed by the same decision: **a ring's inner boundary vertices are the
previous level's outer boundary vertices** — the same entries in the same array,
addressed rather than re-evaluated — and every coarse cell meeting the seam is
split into *three* triangles so that it uses the fine level's midpoint too. It
falls out that at most one edge of any coarse cell can be on the seam: the four
cells at the corners of the hole touch it at a point, so there is no corner case.
The seam is then a set of shared edges with shared endpoints, and there is nothing
left for two evaluations to disagree about.

**A level drops what it cannot resolve.** Keeping a component the grid samples
below its Nyquist does not draw the wave, it draws an artefact with the wave's
amplitude — and the closed form is already in this document: the mesh's own error
at a cell centre is `a (1 - cos(pi / n))` for `n` cells per wavelength, which at
`n = 2` is exactly `a`. **At two cells per wavelength carrying the component is
precisely as wrong as dropping it, and worse to look at, because the artefact has
structure.** Two is therefore the floor and not the answer; the default
`minimumCellsPerWavelength` is **four**, where the residual is 29%. Components are
sorted by wavenumber once per build so each level's set is a prefix, which makes
the sets nested: a component a ring carries is carried by every finer level.

The same `sim::WaveField` still drives every level. A band-limited level is this
spectrum with its unresolvable tail removed, not a second one.

### Cascade — measured cost

At **equal near-field cell size and equal reach** — 2.48 m cells out to 525 m,
128 components (Hs 4 m, Tp 9 s, 16 × 8), one core, `-O2`:

| | cells / levels | vertices | displacement |
|---|---|---|---|
| uniform patch | 424 cells | 180 625 | **55.0 ms** |
| cascade | 2 levels | 79 289 | **25.4 ms** |

Cheaper twice over: three quarters of the outer area at a quarter of the vertex
density, and an outer ring that drops what it cannot carry. Then the interesting
number — the same cascade taken to the horizon:

| | levels | reach | vertices | displacement |
|---|---|---|---|---|
| cascade, matched reach | 2 | 525 m | 79 289 | 25.4 ms |
| cascade, horizon reach | 10 | 134 km | 350 649 | 45.9 ms |

**256 times the reach for 1.8 times the cost.** The uniform grid that reached
134 km at 2.48 m cells would be 1.2 × 10¹⁰ vertices — 33 600 times the cascade,
and about four hours a frame.

`tools/seaway_view` (S-175, Hs 4 m, camera 69 m up at 1280 × 720, GTX 1070 Ti):

| | reach | vertices | displacement | upload | GPU |
|---|---|---|---|---|---|
| uniform patch (before) | 525 m | 180 625 | 55 ms | 1.4 ms | 0.2 ms |
| cascade (now) | 67 km | 316 729 | **46 ms** | 2.5 ms | 0.24 ms |

Nine levels of 212 cells, and what each of them carries is the point:

| level cell size | 2.5 m | 5 m | 10 m | 20 m | 40 m | 79 m and out |
|---|---|---|---|---|---|---|
| components carried | 128 | 128 | 120 | 96 | 16 | **0** |
| vertices | 45 369 | 33 920 | 33 920 | 33 920 | 33 920 | 33 920 each |

The four outermost rings resolve nothing, skip the recurrence entirely and are
dead flat, which is what makes 60 km of sea nearly free — and is correct as well
as cheap, since a 20 m wave at 20 km is a thousandth of a pixel.

**What the cascade did *not* buy.** The near field still costs what it always did:
level 0 and level 1 carry every component and are 79 000 of the 317 000 vertices
but most of the 46 ms. Reaching the horizon is now cheap; resolving the water
under the bow is not, and the FFT above is still the answer to that one. The
vertex buffer also grew 1.75×, and it is re-uploaded whole every frame — 2.5 ms,
which is more than the GPU spends drawing it. The index buffer never changes
between frames and has no business being re-uploaded at all; that is the first
thing to fix if this path stays on the CPU.

### Cascade — limits

- **The rings are squares, so reach is not isotropic.** The corners are √2 further
  out than the sides, and `reach()` reports the side. A camera looking along a
  diagonal has 41% more sea than it needs; the criterion is stated for the worst
  case, which is the axis.
- **The cascade is centred on the ship, not on the camera.** `seaway_view`'s eye
  rides 1.42 Lpp out on the quarter, so level 0 is sized to contain both. A camera
  taken further out than `innerHalfExtent` would find the water directly beneath
  it drawn at level 1's resolution.
- **The sea seen at the horizon is seen edge-on, and a linear facet can then
  present its underside.** Where a far facet tilts away from the eye more steeply
  than the view ray, `gl_FrontFacing` is false, the shared fragment shader flips
  the normal as it must for a hull seen from inside, and the pixel comes out
  black. Measured on `seaway_view`: **9 pixels in 8 frames of 1280 × 720** — about
  one in a million — all within 60 rows below the horizon, which is the far tenth
  of the sea. The same 9 at Hs 4 m and at Hs 5.5 m, so it is not a steepness
  effect waiting to grow. It is a grazing-incidence artefact of a heightfield
  rather than a hole — the elevation channel shows the fragments are drawn, and
  the crack test finds no background there — and the fix belongs with the
  reflection work, where the sea stops going through a hull's BRDF.
- **A far plane has to be moved to match.** Geometry beyond it is clipped, and a
  clip plane makes exactly the edge the cascade exists to remove. `seaway_view`
  now sets far to 1.5 × the reach and near to 0.02 Lpp, which keeps depth
  resolution at the waterline near a millimetre in float32.

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

### Cascade — what the tests are pointed at

A cracked sea and a whole one look identical everywhere except at the crack, and a
crack is a few pixels. So the primary instrument is not a picture at all:

- **An edge census over the whole cascade.** Every directed edge of the triangle
  soup is counted. An interior edge of a watertight surface appears exactly twice
  in opposite directions; an edge on the outer boundary appears once. **A
  T-junction is exactly an interior edge that appears once**, so the crack question
  becomes combinatorial — no camera, no pixels, no tolerance — and the answer is an
  integer: `4n` unmatched edges, all of them on the outer square, and zero
  malformed. Same instrument as the manifold check that caught the hull wound
  inconsistently. It is also the *only* thing that catches a reversed winding in a
  transition cell, which renders as a dark patch and nothing else notices.
  Unstitched, the same census predicts and finds `4n + 6n(levels − 1)`: three
  edges go unmatched at each of the `2n` coarse cells along each seam.
- **The vertex and triangle counts are closed forms** —
  `(n+1)² + (levels−1)(3n²/4 + n)` and `2n² + (levels−1)(3n²/2 + 2n)` — so a
  topology change is an arithmetic disagreement rather than a slightly wrong
  picture, and no vertex may be displaced and then never drawn.
- **The seam is sampled from both sides.** At the midpoint of every coarse cell
  along a seam, the mesh's own elevation a tenth of a millimetre inside and a tenth
  outside must agree to within `2ε` times the field's steepest slope — derived, not
  chosen. The guard is that the same probes on an unstitched cascade report a
  **0.36 m** step, which is the sagitta the stitch is closing: the comparison
  demonstrably can fail. Points standing *exactly* on a level boundary are asserted
  separately, because a level chosen with `<` rather than `<=` sends them one level
  out and the ε probes step right over it.
- **No background pixel below the horizon, over four azimuths.** The elevation
  channel tags every drawn sea fragment with blue 255 and the frame is cleared to
  blue 0, so "is this pixel sea" is an integer question rather than a colour
  comparison. The horizon row is computed from the camera as the vanishing point of
  the sea plane. Four azimuths because a seam is a square and one camera sees two
  of its four sides — a mutation that removed the transition cells on the other two
  went straight through a single view.
- **Two negative controls on that frame**: the uniform patch it replaces, at the
  same cell size and with *more* vertices than the whole cascade, leaves 61 680
  background pixels below the horizon; the same cascade unstitched leaks 480
  through its seams.
- **The near field still passes the load-bearing check**, unchanged: the cascade's
  level 0 renders to 0.6 mm rms and 2.2 mm worst against `WaveField::elevation()`,
  which are the same figures the uniform patch produces, because it is the same
  grid. The comparison and its tolerance construction are now one shared routine so
  there is one place for them to be right.

### What mutation testing changed, second time

Twenty-five deliberate defects in the cascade. Most were caught loudly — the
delegation from a ring to the level inside it made one index tighter trips 11
assertions, a seam midpoint replaced by the coarse corner trips 4, an elevation
scaled by 1.02 trips 4 and is rejected by 374 of 400 samples. **Five escaped
entirely**, and each is now a test that did not exist:

| escaped | now caught by |
|---|---|
| the flat fast path reusing the previous row's accumulator | a level that resolves nothing being **bit-exactly** flat — `z == 0.0f`, normal `(0,0,1)`. It shares the row buffers with every other level, so the stale values are a plausible sea and every other assertion was about levels that do displace |
| `sampleElevation` choosing its level with `<` instead of `<=` | sampling *exactly on* each boundary, including the outer edge. Every other probe was deliberately a hair to one side |
| the index cache key forgetting `stitchSeams` | building two configurations into the *same* surface and holding the result against one built from nothing. Every other test used a fresh object per configuration, so the cache was never asked to notice a change |
| `cellsPerSide` rounding *down* to a multiple of four | a sweep of every resolution from 1 to 40, asserting the rounding is up and that each one is still one watertight sheet. Rounding down is silent and hands back a coarser grid than the resolution criterion asked for |
| the rings' normals replaced by a flat `+z` | a ring's normals checked against the analytic slope of **what that ring carries** — with a guard that it is measurably not the full field's slope. `docs` already records this exact defect escaping once on the uniform patch; the cascade reintroduced the hole by having its own displacement loop |

A sixth was a partial escape worth recording on its own: removing the transition
cells on two sides of every seam was invisible to the rendered frame (the camera
was looking the other way) and obvious to the edge census. That is the general
lesson again — the configuration that makes an assertion exact is often the
configuration that makes it blind — and it is why the frame now asks four times.

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

### Damage — what exists

`engine/gpu/damage.{hpp,cpp}`, `gpu::SceneMesh::appendShip`'s damaged overload,
`engine/gpu/materials/marine.materials`, `tests/test_hull_render.cpp`,
`tools/ram_view`. A damaged ship draws as damaged: plating pushed in where the
zone says it was pushed in, holes where panels tore, and exposed metal round their
edges.

**It is a sibling of `hull.cpp` and it contains no Vulkan.** The whole of it is
mesh arithmetic — refine, displace, cut — so it builds and is checkable on a
machine with no device, which is the same argument that keeps `material.cpp` out
of the Vulkan-gated sources. That is not cosmetic: the instrument that actually
catches cracks here is a combinatorial edge census, and it should not need a GPU
to run.

**It is a rebuild, not a per-frame step.** `buildDamagedHull` produces a body-frame
`DamagedHull` when the damage changes; `appendShip` applies the rigid-body
transform every frame as it already did. The plan above puts the skinning in a
compute shader over a node buffer that never leaves the GPU, and that is still
right for a zone being solved live. It is not right *yet*, because nothing solves a
zone live — `06-roadmap.md` still lists adaptive promotion as the largest thing
outstanding — and a CPU rebuild that costs 2 ms once is not what a frame budget is
spent on. The interface is the place to change it: `HullDamage::displacementAt` is
a pure function of an undeformed position, which is exactly what a vertex or
compute shader would evaluate.

### Joining a deformed patch to an undeformed hull

Two hazards, and they have separate answers. Neither is a tolerance.

- **A crack** — two triangles that shared an edge no longer agreeing where it is.
  Prevented by construction: the displacement is a **pure function of the
  undeformed position**. Two triangles sharing a position — by index, or as the
  coincident-but-distinct vertex records a mirrored hull carries along its seam —
  are handed the same position and get back the same answer. Nothing is displaced
  per triangle, and the deformation adds and removes no vertices.
- **A seam** — a step where the deformed region stops. For a solved zone this is
  prevented by the boundary condition the solve already has: `zone::Edge::Clamped`
  pins the patch perimeter, so the outermost nodes carry displacement **exactly
  zero** and the field reaches zero *inside* its own support instead of being cut
  off at it. `HullDamage::boundaryDisplacement()` reports the largest displacement
  on the patch's boundary loop and the test asserts it is zero while the interior
  maximum is 0.18 m — so this is measured from the solve rather than assumed of it,
  and a patch solved with a free edge would say so.

  **The general rule is that every source must fade to zero at the edge of its own
  support, never be cut off at it**, and that is easier to get wrong than the crack
  is. The membrane tent got it wrong: its falloff ran on the in-plane radius while
  its support was cut by a slab, so the field stepped at the slab's face. Its
  support is a *ball* now, and the instrument that would have caught it — holding
  the field to its own Lipschitz constant — is below.

**The field is the FEM's own, not a fit to it.** `displacementAt` is piecewise
linear over the patch's own elements, each split into two triangles, evaluated by
barycentric weights — which makes it *interpolating*: at a mid-surface node the
weights are exactly (1, 0, 0) and the answer is exactly that node's displacement.
Measured over every node of a solved patch, **2.8 × 10⁻¹⁷ m**, which is an equality
and is asserted as one. Inverse-distance (Shepard) weighting was written first and
rejected: it interpolates only in the limit, it leaves a flat spot at every data
point, and it turns "does the picture agree with the solver" into a tolerance.

Splitting the quad into triangles rather than inverting the bilinear map is also
what keeps it C0 — along a shared element edge both sides are linear in the same
parameter through the same two nodal values, and the diagonal belongs to one
element — so the field is continuous everywhere and identically zero off the patch.
It is zero, not small: the test asserts `1e-12`-exact zero a hull's breadth away,
a compartment down, forty metres forward, and — the check that actually catches
something — at **every node of the patch mirrored twenty metres inboard**, which is
the port shell node for node.

**A patch is a surface embedded in a solid, so a through-thickness guard is not
optional.** A hull wraps; a point on the far side of the ship projects into the
same in-plane cell as a point on the near side. Without the guard, ramming her
starboard side dents her port side identically. Mutation-tested, and the first
version of that test did not catch it — see below.

### Resolution: a dent needs somewhere to happen

The reference ferry's hull is tabulated at 25 stations and 12 waterlines, so her
plating is drawn at about 5 m × 1.5 m. A 0.2 m dent over a 2 m punch on that mesh
moves nothing at all — there is no vertex inside the dent to move. **The plating
has to be refined where the damage is, and refining part of a mesh is the ocean
cascade's problem in a triangle mesh**: a T-junction is a vertex halfway along a
neighbour's edge that the neighbour interpolates straight past.

The answer has the same shape as the cascade's — make the two sides of every edge
decide identically:

- whether an edge splits is a **pure function of its two endpoints**, their length
  against a target size that depends only on position, so the two triangles sharing
  it cannot disagree;
- the midpoint is `(a + b) × 0.5` on the **welded** endpoint positions, which is
  commutative in IEEE arithmetic and therefore bit-identical from both sides, and
  it is interned against the welded endpoint pair so both sides address one vertex
  rather than two that coincide;
- a triangle with one or two split edges is divided into two or three so that it
  **uses the midpoint too**. That is the cascade's transition cell.

The recursion is uniform — every child is examined at depth + 1 whatever template
produced it — so the depth cap cuts both sides of a shared edge at the same depth
and cannot itself open a crack.

The target size is `fineSize + grading × distance` clamped to `coarseSize`, and
**`coarseSize` defaults to larger than any ship**, which is what makes an undamaged
hull cost nothing and a damaged one cost only where it is damaged.

### A tear is a hole, and what that costs

A torn panel is **removed** from the drawn mesh. Not recoloured, not made
transparent, not moved: the triangles are gone, so the pixel shows what is behind
them. The rasterisation already permitted this and it is worth noticing why —
`cullMode` is `NONE` and the fragment shader flips the normal on a back face,
both of which exist in this document because "a hull is seen from inside once it
is cut away", and this is that case arriving. **No pipeline state changes, no
second pass, no sort, no stencil, no `discard`, and no extra vertex attribute.**

What it costs is that the drawn hull stops being closed, so her interior is shaded
and filled where it used to be free. `sim::Ship::hull` is untouched, so buoyancy
still integrates a closed mesh — the hole exists in the render mesh only, the same
separation `breach.hpp` keeps between a hole and the flooding network.

The hole's edge is only as sharp as the refinement, because a triangle is dropped
when its centroid falls inside a torn panel. Measured on a 4 m × 3 m panel at
`fineSize` 0.20 m: **12.146 m² of hole against 12.000 m² of panel**, an error of
1.2%, and the bound is half a triangle along the perimeter. That is why a torn
panel is also a refinement feature — remove it from the feature set and the hole
becomes the shape of the hull's own tessellation.

### Exposed metal is a material, not a shader branch

The plating within `exposedWidth` of a hole takes a material *name*,
`HullPaint::tornEdge`, resolved against the library exactly like the four paint
bands. The shipped set gains `torn_plate_edge`; a mod restates the block or names
a different surface, and nothing recompiles. It is resolved **only on the damaged
path**, so a `.materials` file that predates it still paints an intact ship.

Two things this does not do, both stated rather than pretended away:

- **It is a band, not a blend.** `HullVertex` carries one material index, so two
  surfaces cannot be mixed per pixel. The plastic-strain-driven blend the plan
  above asks for needs a second index and a weight on the vertex and a change to
  `shaders/hull.frag`. The band's area is a closed form — hole perimeter times
  width — and the test asserts it: 3.29 m² against a nominal 3.50 m².
- **A conductor with nothing to mirror is dark.** `torn_plate_edge` is
  `metalness 1.00`, which is right — a fresh fracture is bare steel — and the
  consequence is that with no image-based lighting its radiance is `base × sky`
  and the torn edge reads as a **dark rim** rather than as bright metal. That is
  the same missing term this document already records against `sea_water`'s
  inflated base colour and the sea's absent specular, and it comes right with the
  reflection work. Lowering the metalness to make the picture nicer would be
  putting a wrong number in a data file to compensate for a missing render
  feature, which is what the material set exists not to do.

Also absent: **procedural jagged geometry** on the torn edge. The hole's outline is
the panel's, sampled at the refinement size, which is a straight-ish edge with a
staircase on it rather than torn steel.

### Damage — measured cost

Test ship (796 hull triangles) at 512 × 384 on a GTX 1070 Ti, best of six:

| | triangles | rebuild | scene build | upload | **GPU** |
|---|---|---|---|---|---|
| undamaged | 796 | 0.05 ms | 0.01 ms | 0.05 ms | **0.009 ms** |
| dented and torn | 13 680 | 2.02 ms | 0.24 ms | 0.23 ms | **0.028 ms** |

and `tools/ram_view` on the ferry at 1280 × 720, with the compartments drawn
behind the shell: 27 631 triangles, **0.15 ms** of GPU. Her 1 196 hull triangles
refine to 7 568, of which 3 345 are cut out — 44.8 m² of hole at 4 m/s — in 3.7 ms.
`--frames=N --out=DIR` writes it; the test suite writes `hull_damaged_ship.png`
into `testing::scratchDir()` beside the ship-in-a-sea frame, for the same reason
that one exists.

Three things the table is saying:

- **Undamaged costs nothing, and that is asserted rather than measured.** The
  undamaged path through `buildDamagedHull` returns the ship's own mesh
  **bit-identically** — same vertex bytes, same triangle indices — the damaged
  `appendShip` overload then produces a byte-for-byte identical vertex and index
  buffer, and the rendered frame is byte-for-byte identical. It falls out of the
  construction rather than being special-cased: with no damage the target size is
  never exceeded, so nothing splits, and the displacement is exactly zero. The one
  place that needed care is that **`-0.0 + 0.0` is `+0.0`**, so a zero displacement
  must not go through the addition at all.
- **The cost is the rebuild, not the frame.** 2 ms of refine-displace-cut against
  0.03 ms of GPU. It happens when the damage changes, which at present is once.
- **The GPU is still not where the frame goes.** 17× the triangles for 3× the GPU
  time, and 0.15 ms for a whole damaged ship, her compartments and a sea.

### Damage — what the tests are pointed at

A dented, torn hull looks damaged whatever it is doing — a worse trap than the
sea's, because damage is *supposed* to look irregular, so a wrong answer has cover.
So nothing is eyeballed.

- **A hole is checked by firing a known ray through it**, at a camera looking square
  at the side so every surface in the ray's path is a plane perpendicular to the
  view axis and its clip-space depth is a closed form in its distance alone. Three
  readings of one pixel, and a recoloured panel fails all three:

  | | material id | depth, from the camera alone |
  |---|---|---|
  | intact | her topside paint | 36.5 m |
  | torn, structure behind it | that structure's own material | 40.0 m |
  | torn, nothing behind it | her topside paint again — the **far side, from inside** | 53.5 m |

  The third row is the one that matters: it says the hole is see-through rather
  than see-a-proxy, and the two depth steps are 1700 codes of the 16-bit channel
  apart, which is asserted so the three equalities cannot be satisfied by a frame
  that never changed.
- **Cracks are counted combinatorially**, because a cracked hull and a whole one
  are identical everywhere except at the crack. Every directed edge over welded
  positions: a refined, dented hull with no tear in it has **zero** unmatched and
  zero malformed edges, and is still `sim::isClosedManifold`. With a tear, the
  unmatched count must equal the hole boundary computed **independently from the
  removed triangles** — 234 either way — so an edge the tear exposed and an edge a
  bad refinement opened are told apart rather than lumped together.
- **The rendered consequence is measured too, and it is nearly blind.** A plate on
  her centreplane, in a material nothing outside the ship uses, is invisible through
  intact plating and visible through a crack — the MaterialId channel makes it an
  integer. Four viewpoints around the dent's normal, with the *undamaged* hull
  through the same four cameras as the check that the plate really is hidden.
  **The unstitched control leaks 2 to 8 pixels per view against the census's 840
  unmatched edges.** That is the general lesson from the cascade again, sharper: the
  frame can see the defect, but only just, and only because the plate was put there
  on purpose. Against the ship alone a crack shows the far side of her own hull —
  drawn, because nothing is culled — and the frame sees nothing at all.
- **The refinement criterion is asserted as an inequality over every edge it
  applies to**: no edge inside the graded region may be longer than the target size
  at its own midpoint. 48 420 edges, worst ratio 0.999.
- **A step in the field is a seam, and it is invisible to all of the above.** It is
  a pure function of position, so nothing cracks; the mesh stays a closed manifold,
  so the census is quiet; the cliff is a fold rather than a hole, so the frame sees
  a shaded surface. So the field is held to its own **Lipschitz constant**, which
  for a tent is exact — a cone of depth `d` over a ball of radius `R` cannot change
  faster than `d / R` — over every edge of the drawn mesh and over a sweep across
  the edge of its support. The steepest edge comes out at exactly `d / R`, which is
  the identity a cone satisfies, with a guard that it is not a hundredth of that.
  Applied to the mesh it is also the physical statement: no edge of the plating may
  be stretched by more than the dent's own slope, because a dent stretches steel and
  does not teleport it.
- **A seam in the *input* is not a seam in the output.** A hull carries
  coincident-but-distinct vertex records wherever it has been mirrored or clipped,
  and the two sides of one must split their shared edge into one midpoint rather
  than two. Refining a hull whose lower half has been given its own copies of every
  corner adds **exactly** the same 7 672 vertices as refining the hull that shares
  them, and produces the same triangles.
- **Guards against vacuity throughout**, each because the first version proved
  nothing: that the refinement happened at all (interior edges > 4× the source
  triangle count, deepest split ≥ 4, recursion guard not reached); that the dent
  moved the mesh (0.45 m of a 0.45 m tent); that the patch moved (0.18 m) while its
  boundary did not (0.0); that the frame both bit-identical paths drew is not blank;
  that an intact ship shows **no** exposed-metal pixel anywhere while a torn one
  shows a band of it.

The deformation is held to the solver two ways, and the difference between them is
deliberate. **On the CPU it is an equality** — every mid-surface node reproduced to
2.8e-17 m — because the field interpolates. **Through the frame it is not**, and
cannot be: the depth channel is read at a pixel *centre*, which is not the projected
node, and the sample is the apex of a dent where the plating a few centimetres
either side is measurably shallower. Measured 202 depth codes against 237 predicted
from the camera and the solver, and the assertion is a quarter — enough to reject a
dent drawn at half depth, not enough to pretend the pixel is a point. The exact
statement is the one that never goes through a pixel.

The frame's depth is separately held against **the same triangles rasterised on the
CPU** — `gl_FragCoord.z` is interpolated linearly in screen space, which is exactly
the plane of the projected triangle, so that is an equality to the channel's own
1/65535 and it is an independently written implementation of the depth test.

### Damage — what mutation testing changed

Thirty deliberate defects, in two batches plus a re-check: eighteen across the
whole path, then ten aimed deliberately at the parts the first eighteen had not
touched. Most were caught loudly — a transition template removed trips 5 and 11
assertions, a template that takes the midpoint of the edge that *did not* split
trips 19, a child wound backwards trips 3 including 14 598 malformed edges, the
displacement applied with the wrong sign trips 2, drawing the undeformed mesh trips
2. **Eight escaped on first exposure**, six of them real holes, and every one is
now caught.

**Three were holes in the first batch, and each is now a test that did not exist:**

| escaped | now caught by |
|---|---|
| the through-thickness guard dropped, so ramming her starboard side dents her port side | every node of the patch, mirrored twenty metres inboard along the patch normal. Three hand-picked far points did not catch it, and the reason is worth keeping: one of them *was* the mirror of the impact point and read zero anyway, because the impact landed on a stiffener line whose nodes `Stiffeners::RigidSupport` pins. The probe was in the right place and the node it mirrored was the one node under the punch that could not move. Sweeping every node instead of choosing three is the same lesson as sweeping alignments 1–256 |
| a zero displacement going through the addition anyway | a four-vertex mesh carrying **negative zeros**. `-0.0 + 0.0` is `+0.0`, and no hull in the project has a negative zero in it — a hull tabulated from stations puts its centreline at `+0.0` — so the hazard is real, invisible on every ship here, and needs a mesh built to show it |
| paint bands decided on the deformed hull instead of the undeformed one | a displacement field with a **vertical** component, because a dent normal to a vertical side moves plating in y and never changes its z, so on the flat of a ship's side the two rules agree exactly and nothing can tell them apart. The guard counts how many triangles they disagree about — 957 — so a scene in which they *cannot* disagree fails loudly instead of passing |

**Three were holes in the second batch, and one of them was a defect in the code
rather than only in the tests:**

| escaped | what it was |
|---|---|
| **the tent's falloff replaced by a constant** — plating displaced by the full depth right to the edge of the dent and undisplaced a millimetre further on | The most important escape here, because a step in the field **is** a seam and every instrument was blind to it: the field is a pure function of position so nothing cracks; the mesh stays a closed manifold so the census is quiet; the cliff is a fold rather than a hole so the frame sees a shaded surface. Now caught by holding the field to its own Lipschitz constant — a cone of depth `d` over a ball of radius `R` cannot change faster than `d / R` — over every edge of the drawn mesh and over a sweep across the edge of the support. Looking for that instrument found a **real defect**: the tent's falloff ran on the in-plane radius while its support was cut by a slab, so the field stepped at the slab's face. On the flat of a side the plating never reaches that face, which is why nothing saw it. The falloff now runs on the full 3D distance and the guard is gone, because a ball of radius smaller than a ship's half breadth already cannot reach her other side |
| the exposed-metal band's plane guard dropped, putting an identical band on her port shell | the band's extent checked in the ship's own coordinates against the rectangle the panel was authored as. A single camera cannot make this claim: a band mirrored onto the far side is hidden behind the ship from every viewpoint that can see the hole |
| the hole's area measured on the deformed triangles rather than the undeformed ones | the one scene that is both dented **and** torn. Every other scene has only one of the two, and with no dent the two areas are identical. A hole is a hole in the *plating*, so its area is the undeformed one; the stretch is 0.08 m² of 12.1 m² here and it grows with the dent |

**Two are equivalent mutants, and saying why is the useful part:**

- **The bucket grid's cell size is not a correctness condition.** Fixing it at a
  metre instead of the largest triangle extent is invisible to every assertion,
  which is right: a triangle is registered in *every* bucket its bounding box
  overlaps and the query uses the same grid, so any cell size finds it. The comment
  that used to sit there claimed otherwise and has been corrected — occupancy is
  all the choice buys.
- **The containment epsilon guards against missing a point, not against getting it
  wrong.** Loosening it from 1e-9 to 1e-2 changes nothing, because a mesh node is a
  vertex of *every* triangle that could claim it and the interpolant is linear along
  a shared edge, so whichever triangle answers gives the same value.

Two other results are worth keeping for what they say about the design rather than
about the tests. Removing a patch from the refinement feature set leaves the drawn
dent at **5.4 × 10⁻¹⁸ m** on the ferry's own plating — which is the measurement
behind "a 0.2 m dent on 5 m triangles moves nothing", stated as a number rather
than as an argument. And ignoring the grading, so the target size steps instead of
grading out, trips exactly one assertion and **not** the crack census — correctly,
because the refinement is crack-free whatever the criterion says. The criterion and
the stitching are independent, and the tests treat them that way.

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

### Fire and smoke — what exists

`engine/gpu/smoke.{hpp,cpp}`, `engine/gpu/smoke_gpu.cpp`,
`engine/gpu/shaders/smoke.{vert,frag}`, `tests/test_smoke.cpp`,
`tests/test_smoke_render.cpp`, `tools/smoke_view`.

**The last sentence of the plan above turned out to be the whole of it, and the
first sentence turned out not to apply.** `engine/sim/fire.hpp` is a *two-zone*
model: per compartment it carries a hot upper layer over a cool lower one — two
masses, two internal energies, two species loadings and an interface height. There
are no density or temperature *fields* to raymarch against, so there is no
raymarch. The medium is exactly uniform inside each layer, which makes the
transfer integral closed form:

    L = B (1 − e^(−k d)) + e^(−k d) L_bg

per layer, composed in the order the ray meets them. Along a straight ray through
a horizontally stratified medium that order is unambiguous even where the geometry
chops the path into pieces, because z is monotone in the ray parameter — so two
path lengths and a sign are the whole of it, and the analytic ray–prism
intersection replaces the march entirely. Extinction and emissivity are the same
exponential (Kirchhoff), so **a transparent gas does not glow at any temperature**
and a compartment holding nothing but ambient air composites to the background bit
for bit.

**What the model cannot support, and is therefore not drawn.** Stated here because
a renderer that invents structure the simulation does not have is a lie that will
be believed:

- **No plume.** Heskestad's correlation appears in `fire.hpp` as an entrainment
  *rate* — one scalar, kg/s — and as a mean flame height. The shape that produced
  the correlation was integrated away when it was fitted and is not recoverable.
- **No flame.** `DesignFire` carries a heat release, a base height and a pan
  diameter; there is no flame temperature anywhere in the model. A glowing cylinder
  of Heskestad's mean height would be a scalar wearing a shape.
- **No ceiling jet and no horizontal structure at all.** A zone is well mixed by
  definition.
- **No vertical structure inside a layer**, and the interface is drawn sharp,
  because `fire.hpp` says it is a plane. Softening it would be drawing a mixing
  layer the model does not have — and it would cost the assertion that the
  interface lands where `sim::clipToPixel` puts it.
- **No scattering.** Smoke scatters about as much as it absorbs; treating
  extinction as pure absorption makes a lit layer too dark and stops it glowing
  when a torch is shone into it. In-scattering needs a light list, which is a
  renderer this one is not yet.

**The prism is the model's, not the compartment's.** `Model::attach` sets
`floorArea = gasVolume / height` and says why in as many words: a bounding box
"would over-state the floor area by the turn of the bilge and the interface would
descend too slowly by exactly that ratio". So the drawn volume is that prism —
the compartment's plan aspect ratio scaled to the model's own floor area — and the
drawn gas volume equals `gasVolume`, the drawn hot-layer volume equals
`upperVolume()`, both to machine precision. Measured on the ferry, how box-like
each compartment is (mesh volume as a fraction of its bounding box):

| engine rooms | holds | vehicle deck | accommodation | forepeak | fwd wing tanks |
|---|---|---|---|---|---|
| 99.0% | 96.7–98.7% | 93.4% | 98.8% | 57.8% | 58.0% |

A machinery-space fire is therefore drawn on essentially the right box; a forepeak
fire is drawn on the box the *model* used, which is smaller than the space, and
that is the model's approximation showing rather than the renderer's. The drawn
prism is smaller again by the compartment's permeability — 84.1% of the engine
room's bounding box in plan — because `gasVolume` is the volume of gas and the
machinery displaces the rest.

**Optics from the state, with one knob and one display mapping.** Extinction is a
mass extinction coefficient times a concentration the model states:
`k = K_m · products / V_layer`, with `K_m = 8700 m²/kg` (the SFPE Handbook's
recommended value for the smoke of flaming combustion). `fire.hpp`'s tracer is a
0.05 kg/kg yield, which is a soot-like yield rather than a total-products one, so
that is the right coefficient to apply to it; `SmokeShading::massExtinction` is
the one knob if a ship's fuel says otherwise. Emission is Planck's law integrated
over three contiguous 100 nm bands (R 600–700, G 500–600, B 400–500 nm) — a
spectral binning, not a colourimetric transform, because there is no CIE data in
this repository to do one with. The whole-spectrum integral recovers
Stefan–Boltzmann to 1e-9 and the peak sits on Wien's wavelength to 1e-6, both
asserted.

The single display mapping is `SmokeShading::referenceTemperature`, the gas
temperature whose brightest band reads 1.0. It has to exist because the lit pass's
radiance scale is arbitrary, and it has to be a *temperature* rather than a gain
because the exponential is brutal: **the visible-band radiance of a grey body
rises by a factor of about 400 between 900 K and 1200 K**, so no linear mapping
shows both.

**What a two-zone fire actually looks like, measured.** `smoke_view` on the ferry,
4 MW fast-growth machinery fire in `engine_room_s` (1150 m³, both engine rooms
tracked, watertight door open to port):

| t (s) | Q (MW) | T_upper (K) | interface z (m) | k (1/m) | τ across the room | visibility 3/k (m) |
|---|---|---|---|---|---|---|
| 0   | 0.00 | 288 | 7.00 | 0.00  | 0    | clear |
| 171 | 1.38 | 343 | 3.43 | 2.34  | 60   | 1.28  |
| 343 | 4.00 | 523 | 3.05 | 11.22 | 288  | 0.267 |
| 600 | 4.00 | 531 | 3.13 | 19.90 | 511  | 0.151 |

Three findings, none of them expected:

- **This fire has no fire in it.** A grey layer first puts one byte of red on the
  screen at **834 K**; this one peaks at **531 K** and emits 9.6e-10 of full scale.
  What a 4 MW machinery fire in an engine room looks like, through a two-zone
  model, is *smoke*. The glow is a result and not a setting, and `smoke_view`
  asserts its absence rather than tuning it into view. Raising the design fire to
  10 MW reaches 862 K and the layer does then glow.
- **The layer goes optically black within about a minute and stays there.** τ
  across the room passes 10 before t = 100 s and reaches 511. So the only visible
  information a two-zone model carries is the **interface height** — there is no
  gradient to see, because there is no gradient in the model.
- **The layer does not descend monotonically.** It reaches 2.96 m at about t = 300
  s as the fire stops growing, then *recovers* to 3.13 m as the room settles into
  its vented steady state and the hot layer loses mass through the door. The first
  version of the test asserted a monotone descent and was wrong about the physics,
  not about the code.

**Cost on the target card** (GTX 1070 Ti, 960 × 540, two compartments):
**0.007 ms** for the volumetric pass against 0.056 ms for the lit solid it sits
on. There is nothing to optimise: an analytic two-slab integral is one `exp` per
layer per pixel and the geometry is 36 generated vertices per compartment. Two
passes in one command buffer — the lit solid through `HullRenderer`'s own SPIR-V,
then the medium composited over it with premultiplied alpha, so `exp(−k d)` is the
blender's destination factor rather than a second copy of Beer–Lambert in the
shader.

**How it is checked.** Closed forms on pixels, all agreeing to **one
least-significant bit** — which is the whole budget an 8-bit store allows:

- a uniform slab of known `k` and known thickness transmits `exp(−k t/|dir.x|)`,
  predicted in one line with no geometry code behind it;
- doubling the thickness squares that;
- the medium stops at whatever is solid, at the path length in front of it;
- every pixel against `engine/gpu/smoke.cpp`, which states the same transfer
  integral in double from the formula rather than from the shader — the same
  arrangement `material.hpp` has with `hull.frag`;
- the interface is bracketed by the projections of the interface on the prism's
  near and far faces (2 px apart at 70 m), is exactly the clear colour below the
  first and opaque above the second, and is monotone between them;
- **zero smoke is the unsmoked image, bit for bit** — against `HullRenderer`'s own
  output, which also proves the lit pass here is that pass and not a copy that
  drifted, and with a *drawn* volume of zero extinction as the stronger case.

**One defect this suite would not have found on its own.** The pass first shipped
culling the wrong face. From outside a volume both senses give exactly one fragment
per pixel and exactly the same colour, so every closed-form check passed at 1 LSB;
they differ only with the camera *inside* the medium, where keeping the near faces
draws nothing at all, and that test was off by 152.

**And four more the suite did not find, which mutation testing did.**
`tools/smoke_view/mutate.py` runs 39 single-edit mutants from a copy of the tree
outside the repository, with a per-mutant timeout, a negative control on clean
source and four deliberate equivalent-or-unreachable edits that must survive. The
first pass killed 31 of 35 real mutants with all four controls behaving, and each
of the four survivors was a question no test asked:

- **the layer order.** Every camera in the render suite happened to produce
  ascending rays once the ship's heel was taken out of the body frame, so forcing
  the composite order left the whole suite green. The per-pixel sweep now runs two
  cameras and counts rays of both orders.
- **`expm1` against `exp(x) − 1`.** The comment justifying `expm1` said the
  whole-spectrum integral would fail without it, and that was false — the two
  differ by 1e-16 there, because the region that cancels carries almost none of the
  power. What does see it is the Rayleigh–Jeans *limit*, now asserted as a
  convergence rate rather than at a point. Nothing tests a comment.
- **the prism's aspect ratio.** A square footprint of the same plan area satisfies
  the gas volume, the layer split and the interface height, and draws a 24 m engine
  room as a 15 m square.
- **the cool layer's own volume.** The ferry's lower layer never carries any
  tracer, so sizing its concentration by the whole gas space produces the same
  zero. That question now gets asked on a state set by hand.

After those four checks were added, the survivors and the controls were re-run and
every survivor was killed.

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
