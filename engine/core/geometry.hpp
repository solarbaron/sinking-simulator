// SPDX-License-Identifier: MIT
// Closed-triangle-mesh volume integration with half-space clipping.
//
// Everything hydrostatic in this engine reduces to one question: given a closed
// mesh and a plane, what is the volume and centroid of the part below the plane?
// Buoyancy asks it of the hull against the sea surface. Floodwater asks it of a
// compartment against the internal free surface. Both use integrateBelowPlane().
#pragma once

#include "math.hpp"

#include <cstdint>
#include <vector>

namespace sim {

struct Tri {
    std::uint32_t a, b, c;
};

struct TriMesh {
    std::vector<Vec3> verts;
    std::vector<Tri>  tris;

    void append(const TriMesh& other);
};

// Result of integrating a (sub)region of a closed mesh.
struct VolumeIntegral {
    double volume = 0;    // m^3
    Vec3   centroid{};    // m, valid only when volume > 0
};

// Full closed-mesh volume and centroid. Assumes outward-facing winding.
VolumeIntegral integrate(const TriMesh& mesh);

// Volume and centroid of the region satisfying dot(n, x) <= planeOffset.
//
// No cap geometry is constructed: the divergence-theorem integral is taken about
// a reference point lying *on* the cutting plane, so every cap tetrahedron is
// degenerate and contributes exactly zero. This makes the routine correct for
// meshes of any topology, including ones the plane cuts into several loops.
VolumeIntegral integrateBelowPlane(const TriMesh& mesh, const Vec3& n, double planeOffset);

// Second moment of area of the plane-mesh intersection (the "waterplane"), about
// axes through the intersection's own centroid. i.x is the transverse moment that
// drives free-surface effect; i.y is the longitudinal moment. Returns area in .z.
// Approximated by finite difference of the clipped volume centroid, which is exact
// for prismatic sections and adequate elsewhere.
Vec3 waterplaneMoments(const TriMesh& mesh, const Vec3& n, double planeOffset);

// Solve for the plane offset that makes integrateBelowPlane() return targetVolume --
// i.e. "where does the water surface sit inside this compartment?". Bracketed
// regula falsi; volume is monotone in the offset, so it always converges.
double solvePlaneOffsetForVolume(const TriMesh& mesh, const Vec3& n, double targetVolume,
                                 double loOffset, double hiOffset, int iterations = 40);

// --- Primitive builders -----------------------------------------------------

// Axis-aligned box from min to max corner, outward winding.
TriMesh makeBox(const Vec3& lo, const Vec3& hi);

// A ship hull generated from transverse station offsets, the way real hull forms
// are tabulated. Each station is a half-breadth curve sampled at waterlines; the
// mesh is mirrored about the centreplane and capped at bow, stern and keel.
struct Station {
    double x = 0;                  // longitudinal position, m
    std::vector<double> halfBeam;  // half-breadth at each waterline, m
};
TriMesh makeHullFromStations(const std::vector<Station>& stations,
                             const std::vector<double>& waterlines);

// Transform every vertex of a mesh by R*v + t. Used to bring compartment meshes
// from their authoring frame into the ship body frame.
TriMesh transformed(const TriMesh& mesh, const Mat3& R, const Vec3& t);

}  // namespace sim
