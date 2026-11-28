title: Sample a volume into hexahedra
category: Assembly and analysis

The public `create` binding takes sample-point counts along U, V, and W. Studio specifies the resolution in cells. An axis with \(m\) cells needs \(m+1\) sample points. For \(4\times3\times2\) cells, request \(5\times4\times3\) points.

Start with the [cube example](https://github.com/gha3mi/forcad/blob/dev/example/cmp_volume.f90). Add these declarations beside its `Xc` declaration.

```fortran
real(rk), allocatable :: points(:, :)
integer, allocatable :: cells(:, :)
```

After `shape%set(...)`, replace the volume measurement and its print statement with the following calls.

```fortran
call shape%create(res1=5, res2=4, res3=3)
if (.not. shape%err%ok) error stop "Could not sample the volume"
points = shape%get_Xg()
cells = shape%get_elem_Xg_vis()
print '(4(i0,1x))', size(points,1), size(points,2), size(cells,1), size(cells,2)
```

The program prints `60 3 24 8`. It produced 60 points with three coordinates each and 24 cells with eight node indices each. Each index identifies a row in `points`. To inspect the grid outside your program, [write it to VTK](export-volume-vtk.html).

**Studio equivalent.** [Mesh a single volume](../how-to/mesh-a-volume.html) with 4, 3 and 2 cells, then check the node and cell counts in the status bar.

**Multipatch meshes.** Each member has its own sampled grid. For a single, conforming mesh across joined faces, [generate it in Studio](../how-to/mesh-a-multipatch.html). To inspect the complete reconstruction in Fortran, [export the project as a program](../how-to/export-fortran.html). The exported mesh recipe includes the step that assigns common node numbers across the joined faces.
