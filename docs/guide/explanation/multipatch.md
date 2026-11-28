title: Multipatch interfaces and conforming meshes

## What the model stores

Each member of a multipatch model retains its control points, weights, and knot vectors. A connection records the two boundaries and how their parameters correspond. For example, a volume face can meet another face with the two tangential directions exchanged or reversed.

With positional continuity (\(C^0\)), corresponding points on the two boundaries coincide. A \(G^1\) join also relates first derivatives through a change of parameters across the interface. For two regular surface patches, this gives a common tangent plane along the seam. Moving a control or changing its weight can therefore require changes to neighboring controls. Studio previews those dependencies in the [shared control editor](../how-to/edit-shared-control.html).

## What the mesh shares

For a conforming hexahedral mesh, the sampled nodes on joined faces must match after applying the interface mapping. If the first face’s U direction maps to the second face’s V direction, eight cells along U on the first face require eight along V on the second.

The public multipatch sampling call constructs a grid for each patch. When [generating a conforming mesh](../how-to/mesh-a-multipatch.html), Studio assigns the matching nodes a common global ID. The [exported Fortran program](../how-to/export-fortran.html) includes this numbering step.

## What an analysis shares

In a scalar IGA problem, degrees of freedom are coefficients of the solution field. On a matching clamped face with compatible trace weights, `cmp_elem(shared=.true.)` gives the two elements common IDs for their face coefficients. The geometry and its visualization can still keep separate arrays for each member. Changing the visualization sampling density leaves those solution IDs unchanged.

The shared coefficients impose equality of the field traces. Flux balance has a different sign convention because the two outward normals point in opposite directions. The [two-volume Poisson guide](../iga/solve-two-volumes-poisson.html) checks each flux against a known solution and integrates it over the face. For a geometric design variable, [differentiate the volumes on either side of a shared face](../fortran/shared-face-volume-sensitivity.html) to see how the two patch contributions combine.
