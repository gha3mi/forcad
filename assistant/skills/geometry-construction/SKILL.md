---
name: geometry-construction
description: Plan and build exact ForCAD curves, surfaces, and volumes from dimensions, sketches, or reference attachments.
---

# Construct geometry

1. Identify the requested dimensionality, feature relationships, and tolerances. Read the model units. Treat a photo or drawing as evidence for visible form and marked dimensions. Flag hidden geometry or scale that must be inferred.
2. Use `build_model` for a sequence of native primitives and construction operations that can refer to earlier steps in the same call. Use `construct` for extrude, revolve, sweep, or loft from existing patch IDs. Consider `build_curved_volume` when the requested solid fits its axial NURBS multipatch representation. For general multipatch designs, plan native pieces and their interfaces.
3. Explain which steps establish the exact shape and where the inputs include estimates. After building, inspect the resulting patch or assembly, measure the relevant dimensions, and check diagnostics. Confirm that a requested solid is represented by volume patches.
