title: Export a sampled volume to VTK
category: File exchange

Write a 12 × 8 × 4 mm block as an unstructured grid in physical space. The three sample counts control the exported resolution. The coordinates in the VTK file use the same length unit as the input.

Use the [fpm project from the first curve tutorial](../tutorials/first-curve.html#add-the-dependency) and replace `app/main.f90` with this program.

```fortran
program export_sampled_volume
    use forcad, only: nurbs_volume, rk
    implicit none
    type(nurbs_volume) :: block

    call block%set_hexahedron(L=[12.0_rk, 8.0_rk, 4.0_rk], nc=[2, 2, 2])
    if (.not. block%err%ok) error stop "Could not construct the block"
    call block%create(res1=4, res2=3, res3=2)
    if (.not. block%err%ok) error stop "Could not sample the block"
    call block%export_Xg("block.vtk", encoding="ASCII")
    if (.not. block%err%ok) error stop "Could not export the sampled block"
end program export_sampled_volume
```

Run `fpm run --profile release` from the project directory. Open `block.vtk` in a text editor and check these records.

| Record | What to check |
| --- | --- |
| `POINTS 24 double` | The \(4\times3\times2\) samples span X from 0 to 12, Y from 0 to 8, and Z from 0 to 4. |
| `CELLS 6 54` | The grid contains \(3\times2\times1=6\) cells. Each row has an `8` followed by eight point indices. |
| `CELL_TYPES 6` | The six following values are all `12`, the VTK identifier for a linear hexahedron. |

Choose `encoding="BINARY"` for binary output once you have checked the ASCII example. For the same export in Studio, use [Export sampled geometry as VTK](../how-to/export-geometry-vtk.html). To export the nodes and cells of an existing Studio mesh, [export its mesh row](../how-to/export-mesh-vtk.html).
