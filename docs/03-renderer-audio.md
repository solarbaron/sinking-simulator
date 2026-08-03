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
