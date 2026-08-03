// SPDX-License-Identifier: MIT
//
// Validation of the camera matrices.
//
// These exist so that every rendering assertion downstream can be stated as a
// closed form -- "this world point lands on that pixel" -- rather than as
// "the picture looks about right". Which means they have to be right first.
//
// The failures worth guarding against are the ones that produce an image that
// looks almost correct: a Y axis flipped for OpenGL instead of Vulkan gives a
// vertically mirrored render, and an OpenGL [-1,1] depth range on Vulkan
// hardware silently uses half the depth buffer.
#include "engine/core/math.hpp"
#include "harness.hpp"

#include <cstdio>

using sim::Mat4;
using sim::Vec3;
using testing::expectNear;
using testing::expectTrue;

namespace {

void testMatrixAlgebra() {
    const Mat4 identity = Mat4::identity();
    double out[4];
    identity.transform({3, -4, 5}, out);
    expectNear("identity leaves x", out[0], 3.0, 1e-12);
    expectNear("identity leaves y", out[1], -4.0, 1e-12);
    expectNear("identity leaves z", out[2], 5.0, 1e-12);
    expectNear("identity gives w = 1", out[3], 1.0, 1e-12);

    // Multiplication must be associative and identity-preserving; a transposed
    // operator() would still pass a symmetric test, so use an asymmetric matrix.
    Mat4 a = Mat4::identity();
    a(0, 3) = 10.0;  // translation in x
    Mat4 b = Mat4::identity();
    b(1, 3) = -7.0;  // translation in y
    const Mat4 ab = a * b;
    ab.transform({0, 0, 0}, out);
    expectNear("composed translation applies x", out[0], 10.0, 1e-12);
    expectNear("composed translation applies y", out[1], -7.0, 1e-12);
    expectTrue("identity is a left unit", (Mat4::identity() * a).m == a.m);
    expectTrue("identity is a right unit", (a * Mat4::identity()).m == a.m);
}

// Vulkan clip space: z in [0, 1], not OpenGL's [-1, 1].
void testPerspectiveUsesVulkanDepthRange() {
    const double nearPlane = 0.5, farPlane = 100.0;
    const Mat4 projection = sim::perspective(60.0 * sim::kDegToRad, 16.0 / 9.0, nearPlane, farPlane);

    double clip[4];
    projection.transform({0, 0, -nearPlane}, clip);
    expectNear("a point on the near plane maps to depth 0", clip[2] / clip[3], 0.0, 1e-9);

    projection.transform({0, 0, -farPlane}, clip);
    expectNear("a point on the far plane maps to depth 1", clip[2] / clip[3], 1.0, 1e-9);

    // Halfway in *depth* is not halfway in distance -- perspective depth is
    // nonlinear, and a linear result would mean an orthographic matrix slipped in.
    projection.transform({0, 0, -(nearPlane + farPlane) * 0.5}, clip);
    const double midDepth = clip[2] / clip[3];
    expectTrue("perspective depth is nonlinear", midDepth > 0.9 && midDepth < 1.0);
}

// The y flip is the difference between a correct render and a vertically
// mirrored one, which is easy to miss by eye on a symmetric scene.
void testPerspectiveFlipsYForVulkan() {
    const Mat4 projection = sim::perspective(90.0 * sim::kDegToRad, 1.0, 1.0, 100.0);
    double clip[4];

    // A point above the eye must land in the *upper* half of the image, which in
    // Vulkan means negative clip y.
    projection.transform({0, 1, -1}, clip);
    expectTrue("a point above the axis has negative clip y", clip[1] / clip[3] < 0.0);

    projection.transform({0, -1, -1}, clip);
    expectTrue("a point below the axis has positive clip y", clip[1] / clip[3] > 0.0);

    // At 90 degrees vertical field of view and unit aspect, the frustum edge at
    // one unit of depth sits at exactly one unit of height.
    projection.transform({0, 1, -1}, clip);
    expectNear("the vertical frustum edge maps to the clip boundary", clip[1] / clip[3], -1.0,
               1e-9);
    projection.transform({1, 0, -1}, clip);
    expectNear("the horizontal frustum edge maps to the clip boundary", clip[0] / clip[3], 1.0,
               1e-9);
}

void testLookAtPlacesTheCamera() {
    // Camera at +Z looking back at the origin, z-up world.
    const Mat4 view = sim::lookAt({0, 0, 10}, {0, 0, 0}, {0, 1, 0});
    double out[4];

    view.transform({0, 0, 10}, out);
    expectNear("the eye maps to the view-space origin (x)", out[0], 0.0, 1e-12);
    expectNear("the eye maps to the view-space origin (y)", out[1], 0.0, 1e-12);
    expectNear("the eye maps to the view-space origin (z)", out[2], 0.0, 1e-12);

    // The target is 10 units in front, which is -Z in view space.
    view.transform({0, 0, 0}, out);
    expectNear("the target sits in front of the camera", out[2], -10.0, 1e-12);

    // A point to the world's +X must appear to the camera's right when looking
    // down -Z, i.e. view-space +X.
    view.transform({3, 0, 0}, out);
    expectNear("world +x maps to view +x", out[0], 3.0, 1e-12);

    // lookAt must not depend on the length of the up vector.
    const Mat4 longUp = sim::lookAt({0, 0, 10}, {0, 0, 0}, {0, 50, 0});
    expectTrue("up vector length does not matter", longUp.m == view.m);
}

// The closed form every rendering test will lean on: a known world point must
// land on a known pixel.
void testWorldPointLandsOnKnownPixel() {
    constexpr double width = 800, height = 600;
    const Mat4 view = sim::lookAt({0, 0, 10}, {0, 0, 0}, {0, 1, 0});
    const Mat4 projection =
        sim::perspective(90.0 * sim::kDegToRad, width / height, 0.1, 1000.0);
    const Mat4 viewProjection = projection * view;

    double clip[4], x = 0, y = 0;

    // The point the camera is aimed at must land dead centre.
    viewProjection.transform({0, 0, 0}, clip);
    expectTrue("the target projects in front of the eye", sim::clipToPixel(clip, width, height, x, y));
    expectNear("the target lands on the horizontal centre", x, width * 0.5, 1e-6);
    expectNear("the target lands on the vertical centre", y, height * 0.5, 1e-6);

    // Ten units away with a 90 degree vertical field of view, the top edge of the
    // frustum is ten units up. So (0, 10, 0) must land exactly on the top edge.
    viewProjection.transform({0, 10, 0}, clip);
    expectTrue("the frustum top projects", sim::clipToPixel(clip, width, height, x, y));
    expectNear("the frustum top edge maps to pixel row zero", y, 0.0, 1e-6);
    expectNear("the frustum top edge stays horizontally centred", x, width * 0.5, 1e-6);

    // And a point behind the camera must be reported rather than silently
    // wrapping round to a plausible-looking pixel.
    viewProjection.transform({0, 0, 20}, clip);
    expectTrue("a point behind the eye is rejected", !sim::clipToPixel(clip, width, height, x, y));
}

}  // namespace

void runCameraTests() {
    std::printf("\n--- camera matrices ---\n");
    testMatrixAlgebra();
    testPerspectiveUsesVulkanDepthRange();
    testPerspectiveFlipsYForVulkan();
    testLookAtPlacesTheCamera();
    testWorldPointLandsOnKnownPixel();
}
