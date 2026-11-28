title: Project and exchange commands
category: Projects and exchange

## Save or combine projects

[Save the active project](save-project.html), [import an FCAD project into it](import-fcad.html), or [copy selected objects between tabs](copy-between-projects.html).

## Choose an export

<figure class="doc-figure" markdown="1">

<a href="../../media/project-exchange-ribbon.webp">
  <picture>
    <source media="(max-width: 32rem)" srcset="../../media/project-exchange-ribbon-detail.webp" width="288" height="105">
    <img src="../../media/project-exchange-ribbon.webp" width="568" height="105" alt="Exchange group on the Project ribbon, with Export VTK and Export IGES beside Export Fortran">
  </picture>
</a>

<figcaption>On the <strong>Project</strong> tab, geometry exports use your selection. <strong>Export Fortran</strong> uses the entire active project.</figcaption>
</figure>

- **Selected geometry** — [Export sampled geometry as VTK](export-geometry-vtk.html) or [export a curve or surface as IGES](export-iges.html).
- **Whole project** — [Export a Fortran program](export-fortran.html) to rebuild the model.
- **Generated mesh** — Select its row under **Scene → Meshes**, then use [Mesh → Export mesh (VTK)](export-mesh-vtk.html) to write its nodes and cells.

## Import curves or surfaces from IGES

[Import IGES geometry](import-iges.html) covers supported spline entities and source units. See [File formats](../reference/file-formats.html) for the full list of formats.
