# 02 — Simulation

The physics plan, in dependency order. Section 1 is implemented; the rest is the
build.

---

## 1. Hydrostatics, flooding and air — **implemented**

See `README.md` for behaviour and `engine/sim/ship.cpp` for the code. Summary of
what exists and what it still needs.

**Done:** closed-mesh volume/centroid integration under an arbitrary plane;
internal free surfaces re-levelled per tick; floodwater as real mass at its real
centroid (free-surface effect emergent, validated against ρ·i/Δ); compressible
trapped air per compartment with isothermal state and pressure-limited transfer;
a single two-phase orifice network covering breaches, doors, hatches, vents,
pipes and cross-flood ducts; bilge and ballast pumps with head-dependent output;
6-DOF rigid body with measured hydrostatic stiffness, added mass, quadratic drag
and modal damping; damaged GZ curves computed by forced-heel sweep with the
floodwater free to re-level.

**Next, in order:**

1. ~~**Compartment geometry from CSG rather than authored boxes.**~~ **Done.**
   `clipByPlane` / `clipToBox` carve each compartment out of the hull interior:
   Sutherland-Hodgman per triangle, cut edges welded on a spatial hash, chained
   into loops and capped by ear clipping. Validated by cutting a hull into a grid
   and checking the cell volumes sum to the whole. Still to do: exact predicates
   (Shewchuk) for robustness on degenerate input, nested cap loops (a plane
   through a hollow mast), and non-planar bulkheads.
2. **Permeability from contents, not a constant.** Currently a scalar per
   compartment. Should be derived from the actual cargo, machinery and outfit
   volumes placed in the space, and should change as cargo shifts or burns.
3. **Sloshing.** The quasi-static horizontal free surface is right for slow
   flooding and wrong for a half-full tank in a seaway. Escalation path in §6:
   the lumped model stays until a compartment's excitation exceeds a threshold,
   then it is promoted to a particle solver and demoted when it settles.
4. **Downflooding angle bookkeeping.** Track the heel at which each opening
   immerses, and surface it — this is the number that decides survival and it is
   currently only implicit.
5. **Progressive flooding through the structure.** Requires §3: once elements
   tear, new orifices appear automatically.

---

## 2. Seakeeping and hydrodynamics

The ship currently floats in still water. This adds the sea.

### Wave field

Directional spectrum (JONSWAP for fetch-limited, Pierson–Moskowitz for fully
developed, plus swell trains superposed), sampled as several hundred components
and evaluated by FFT on a tiled grid for rendering. For physics, the analytic
sum is evaluated directly at the points that matter — hull panels, openings, the
free surface at each breach — because interpolating a rendering grid loses exactly
the detail that decides whether a hatch immerses.

Nonlinear extensions: second-order Stokes correction for steep seas, and a
breaking criterion driving spray/foam emission.

### Radiation and diffraction

The honest way to do wave loads on a large ship is potential flow. Plan:

- **Offline:** solve the boundary-element (panel) problem over the hull for a grid
  of frequencies and headings, producing added mass A(ω), radiation damping B(ω),
  and diffraction/Froude–Krylov exciting forces. Use NEMOH or Capytaine as the
  reference solver; ship the coefficient tables as ship assets.
- **Runtime:** Cummins impulse-response formulation. The radiation force has
  memory — it depends on the history of motion — so the convolution integral is
  approximated by a fitted state-space model (4–8 states per DOF, identified
  offline). This is the standard marine-simulation approach and it is cheap at
  runtime.
- **Nonlinear Froude–Krylov and restoring**: integrated over the *instantaneous*
  wetted surface each tick rather than the mean position. The engine already has
  exactly the routine this needs. This is what captures a ship's behaviour in
  large waves, where linear theory quietly stops being true.

### Viscous roll damping

Potential flow gives no roll damping at all — it is entirely viscous and it is the
difference between a ship that rolls 8° and one that rolls 35°. Ikeda's empirical
method (friction, eddy, lift, bilge keel, and appendage components) as the
baseline, with per-ship tuning against roll decay tests.

### The hard cases

- **Green water on deck**: waves over the bulwark, water loose on deck. Handled by
  the same escalation as sloshing (§6): the deck becomes a shallow-water solver
  when wet.
- **Slamming**: bow flare and bottom slam impacts, Wagner/von Kármán added-mass
  impact theory to get the pressure pulse, delivered to the FEM as a load.
- **Parametric roll**: emerges naturally from nonlinear restoring in head seas
  once the above is in place. A good validation target precisely because it is
  hard to get by accident.
- **Broaching and surf-riding** in following seas.

---

## 3. Structure — adaptive tetrahedral FEM

The chosen fidelity target is full 3D tetrahedral finite elements. The physical
obstacle is arithmetic: a 300 m hull with 20 mm plating needs ~10 mm elements to
resolve bending through the plate thickness, which is on the order of 10¹¹
elements for the whole ship. That is not a hardware problem to be solved by
waiting.

The resolution is **spatial adaptivity with substructuring**, which is how crash
and blast FEA has always handled the same conflict:

### Three-tier model

**Tier 0 — global girder, always on.** The hull as a Timoshenko beam with
section properties derived from the real structural arrangement. ~200 DOF.
Gives hull girder bending, hogging/sagging in waves, shear, torsion. Runs at
100 Hz for free.

**Tier 1 — reduced 3D, always on.** The full 3D shell/solid model of the entire
ship, condensed offline by **Craig–Bampton component mode synthesis**: retain the
interface DOF plus the lowest few hundred fixed-interface normal modes, discard
the rest. Reproduces the linear-elastic response of the full model to within a
few percent for the frequency range that matters, at 10⁻⁵ of the cost. Everything
away from damage stays here forever.

**Tier 2 — full nonlinear tet FEM, adaptive.** Around an impact, a fire, or a
growing crack, a region is *promoted*: the reduced model is replaced by genuine
3D tetrahedra at full resolution, coupled to the surrounding Tier-1 model through
the retained interface DOF. Inside that region the physics is uncompromised —
real stress tensors, real plasticity, real fracture. Budget 10⁵–10⁶ elements
across all active zones.

This is not a compromise on the answer where the answer matters. A collision
loads a 20 m stretch of side shell; the response of the bow 200 m away is linear
and the reduced model gets it right. The nonlinearity is local, so the expensive
model is local.

### Constitutive model

- **Elasticity**: co-rotational formulation to handle large rotations with small
  strains; St. Venant–Kirchhoff or Neo-Hookean where strains get large.
- **Plasticity**: J2 (von Mises) flow with combined isotropic and kinematic
  hardening; radial return mapping. Rate dependence via Johnson–Cook, because
  steel is markedly stronger under impact loading and ignoring that overestimates
  damage.
- **Damage and fracture**: ductile damage accumulation (Johnson–Cook damage or
  Gurson–Tvergaard–Needleman for void growth, which is more physical for the
  thick sections that matter). Failure surface reached → element splits along the
  maximum principal stress plane. Element deletion is the cheap fallback but
  loses mass and looks wrong; the plan is mesh splitting with remeshing of the
  affected neighbourhood.
- **Buckling**: falls out of the geometrically nonlinear formulation provided the
  mesh resolves the plate panels and the initial imperfections are seeded.
  Stiffened panel collapse is the failure mode that actually breaks ships' backs
  and it must not be scripted.
- **Welds and joints**: modelled as cohesive-zone interface elements with their
  own (lower) strength. Structures fail at connections far more often than in the
  middle of a plate.

### Element technology — revised after measurement

The original plan was uniform linear tetrahedra throughout Tier 2. **The spike in
`07-fem-spike-findings.md` ruled that out**, and the reason is worth stating
because it is not the obvious one.

Linear tets lock in bending: measured error against beam theory is 63% at two
elements through the thickness, 32% at four, 11% at eight. Ship plating is thin,
so avoiding that error needs many elements through 20 mm of steel — and the
explicit stability limit is set by the *smallest* element dimension, so those
same elements collapse the timestep. Cost scales as h⁻⁴. The two constraints
close on each other and leave no workable resolution.

So Tier 2 is mixed:

- **Solid-shell / assumed-strain (EAS, ANS) elements for plating.** One element
  through the thickness, no locking, and the timestep is governed by in-plane
  size instead of plate thickness — worth 5–10× on the step alone, on top of the
  element count reduction.
- **Tetrahedra where the geometry really is three-dimensional**: castings, engine
  seats, thick brackets, and the crush zone once plating has folded and shell
  kinematics no longer apply.
- **Promotion from shell to tet** as an element crumples past that point.

### Solver

Explicit central-difference time integration inside Tier 2 (standard for
impact/fracture, no global stiffness matrix, trivially parallel). Tier 1 is
implicit and cheap. Lumped mass matrix.

GPU: the element loop is a good compute-shader workload — gather nodal state,
compute deformation gradient, stress, internal force, write per-element forces,
then gather into nodes through a CSR adjacency. **Measured on a GTX 1070 Ti:
450–670 M element-updates/s**, roughly 100× a single CPU core and 4× the whole
24-thread CPU. Nodal forces are gathered rather than scattered, which avoids
float atomics entirely and fixes the accumulation order.

### Slow damage

Separate from the fast solver, on a very long timescale:

- **Fatigue**: rainflow counting on the Tier-0/Tier-1 stress history, Miner's rule
  against S-N curves per detail class. A ship that has worked hard in heavy
  weather for years has cracks where the hot spots are.
- **Corrosion**: thickness diminution by zone, driven by coating condition,
  ballast/cargo/atmosphere exposure, and time. Directly reduces section modulus
  and therefore the loads at which everything above triggers.

### Material database

Per material: density, E, ν, yield and ultimate strength, hardening curve,
fracture strain vs triaxiality, Johnson–Cook rate and thermal coefficients,
thermal conductivity, specific heat, expansion, melting point, and — critically —
**temperature-dependent strength reduction curves** (Eurocode 3 for structural
steel, which loses roughly half its yield at 550 °C and nearly all of it at
800 °C). Coverage: mild steel, higher-tensile grades AH/DH/EH 32/36/40, stainless,
aluminium 5083/5383/6082 (which loses strength at *200* °C — the reason aluminium
superstructures are a fire problem), GRP and sandwich laminates, timber, ferro-
cement, and HY-80/100 for naval hulls.

---

## 4. Fire, heat and gas

Fire matters because of what it does to the structure and the atmosphere, not
because of the flames.

### Combustion

Two-tier, like the FEM:

- **Baseline: multi-zone model.** Each compartment carries an upper hot layer and
  a lower cool layer with an interface height — the classical CFAST formulation.
  Cheap, well validated, correct for smoke filling and layer descent, which is
  what kills people.
- **Promoted: LES combustion** in the compartment of interest. Eddy-dissipation or
  a flamelet model on a coarse (10–20 cm) grid. This is FDS-class physics at
  game resolution, used only where the player is.

Fuel is real: each compartment has an inventory of combustibles with heat of
combustion, ignition temperature, pyrolysis rate and oxygen demand. Fire goes out
when it runs out of either, and an under-ventilated fire produces carbon monoxide
instead of carbon dioxide — which is the actual hazard.

### Gas transport

Species-resolved (O₂, CO₂, CO, unburnt fuel, water vapour, soot) advected through
the same opening network the water uses. A door is a door. Buoyancy-driven flow
through vertical openings uses the standard two-way orifice formulation (hot gas
out the top, cool air in the bottom, with a neutral plane in between).

### Heat

- Conduction through bulkheads and decks by an implicit FEM thermal solve on the
  structural mesh — the same mesh, so the temperature field maps directly onto the
  strength reduction in §3.
- Radiation between hot surfaces and to flame volumes (view factors precomputed
  per compartment).
- Convection to gas layers.

The coupling that makes this worth the effort: **a fire in a machinery space heats
a bulkhead, the bulkhead loses strength, the bulkhead fails under hydrostatic
load from the flooded space next door, and the flooding spreads.** Every step of
that is modelled by a different subsystem and none of them know about the others.

### Suppression

Water spray and deluge (droplet evaporation cooling gas and wetting fuel), CO₂ and
inert gas total flooding (oxygen displacement, with the compartment sealing
requirement that makes it fail if a door is open), high-expansion foam, dry
powder, and boundary cooling. All of them add water or gas mass to compartments,
which is to say all of them affect stability. Firefighting has sunk ships.

---

## 5. Aerodynamics and wind

- **Wind field**: mean profile with a logarithmic boundary layer over the sea
  surface, plus a gust spectrum (NPD or Harris), plus local disturbance around
  the superstructure.
- **Windage**: offline panel-method or RANS sweep over heading angles produces
  force and moment coefficient tables per ship; runtime interpolates. Heeling
  moment from wind is a required input to the IMO weather criterion and to any
  honest capsize model.
- **Superstructure flow**: funnel exhaust dispersion, helicopter deck turbulence,
  and the recirculation zones that decide where smoke goes.
- **Sailing vessels**: for anything with sails, a lifting-line or vortex-lattice
  solver over the sail plan with real angles of attack, sail trim, twist,
  reefing, and the coupled heel/leeway equilibrium. A square-rigger and a modern
  sloop are different aerodynamic problems and both should work.
- **Ventilation systems**: fans, ducts and dampers as another set of edges in the
  gas network, which is how smoke gets somewhere it should not be.

---

## 6. Free-surface fluids

The lumped compartment model is right for slow flooding and wrong for violent
water. The plan is not to replace it but to **escalate**.

Each compartment carries a state: `Quiescent` (lumped, analytic free surface) or
`Dynamic` (particles). Promotion triggers on lateral acceleration, fill fraction
in the sloshing-sensitive band, a nearby impact, or player proximity. Demotion
happens when kinetic energy drops below a threshold for a sustained period, at
which point the particle mass and momentum are integrated back into a lumped
level. Mass is conserved exactly across both transitions by construction.

Solver choice for the dynamic state: **FLIP/APIC on a sparse grid** for interior
water, because it handles the free surface, the pressure projection and the
coupling to compartment boundaries well, and because sparse grids match the
geometry (water is in a few rooms, not everywhere). SPH is the alternative and is
better for spray and jets; the current expectation is FLIP for volumes and a
separate SPH/particle system for spray, jets from breaches, and rain.

Uses: sloshing in partly filled tanks (a genuine stability hazard, and the design
driver for LNG carriers), green water on deck, water moving between rooms as the
ship rolls, the jet from a hull breach, bilge water, and the visual of a
compartment filling around you.

---

## 7. Propulsion, manoeuvring and machinery

- **Propellers**: open-water curves (K_T, K_Q vs J) per propeller, wake fraction
  and thrust deduction from the hull, four-quadrant data so astern and crash-stop
  work, cavitation inception and its thrust breakdown, and ventilation when the
  propeller emerges in a seaway.
- **Manoeuvring**: MMG-style modular model — hull, propeller, rudder terms
  separately identified — rather than a monolithic derivative set, so damage to
  one part degrades one term. Shallow water effects, bank effect, ship-to-ship
  interaction, and current.
- **Steering**: rudder with stall, hydraulic steering gear with real rates and
  failure modes, azimuth thrusters, Voith-Schneider, waterjets, bow and stern
  thrusters with their speed-dependent effectiveness loss.
- **Prime movers**: diesel engine model with turbocharger lag, fuel rack, governor,
  cooling and lubrication circuits, and the failure modes that follow when those
  circuits are damaged or flooded. Steam plant (boiler, turbine, condenser) for
  older and naval vessels. Gas turbines. Diesel-electric and hybrid drivetrains
  with a real electrical load-flow model.
- **Ship systems as networks**: electrical, hydraulic, fuel, fresh water, bilge,
  ballast, fire main, compressed air — each a graph of sources, sinks, valves and
  pipes. Damage severs edges. Loss of the fire main because a pump room flooded
  is a real and common cascade, and it should be modelled as a cascade rather
  than as a scripted consequence.

---

## 8. Cargo, ice and the rest

- **Cargo**: bulk cargo as a granular material that can shift and liquefy (the
  mechanism behind a long list of bulk carrier losses); containers with real
  lashing loads that part in heavy rolling; RoRo vehicles with tyre friction and
  lashings; liquid cargo with sloshing; heavy lift with the crane-induced
  stability problem.
- **Ice**: hull-ice interaction loads for icebreaking, ice-strengthened scantlings,
  and **topside icing** — spray freezing on the superstructure adds high weight
  and has capsized fishing vessels within hours.
- **Grounding**: seabed contact with soil mechanics for the reaction, hull raking
  damage as a moving FEM load, and the stability problem of being partly supported
  by the ground.
- **Collision**: two deformable ships, both FEM-active in the contact zone.
- **Submarines**: pressure hull with depth-dependent loading and collapse depth,
  main ballast and trim tanks, the fact that submerged stability is a different
  problem (no waterplane, so BM = 0 and only BG matters).
- **Crew and passengers**: damage control parties that take real time to reach a
  door, firefighting teams with air supply limits, and passenger evacuation on a
  social-force model over a ship that is listing — walking speed as a function of
  deck angle is a well-measured quantity and it collapses beyond about 20°.

---

## Which of these are load-bearing

If time runs short, the ordering that preserves the most of the concept:

1. Flooding + air + stability **(done)**
2. Waves and seakeeping — without a seaway, nothing else has excitation
3. Structural FEM — the deformation and tearing the concept is named for
4. Fire and thermal — the second casualty type, and the best coupling story
5. Fluids — mostly a fidelity and presentation upgrade over the lumped model
6. Everything else — breadth, and breadth is what makes it a *ship* simulator
   rather than a *sinking* simulator
