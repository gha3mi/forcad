<div class="doc-intro" id="model-and-mesh-spline-geometry">
<div class="intro-copy">
<h1 id="nurbs-geometry-and-meshing"><span>NURBS</span> geometry and meshing</h1>
<p>Model curves, surfaces, and volumes in ForCAD Studio. Use the same geometry in Fortran to fit measured data, mesh joined patches, or assemble an isogeometric analysis (IGA) problem.</p>
<div class="intro-actions">
<a href="page/tutorials/first-curve.html">First curve in Fortran <span aria-hidden="true">→</span></a>
<a href="page/tutorials/explore-studio.html">Explore Chess in Studio <span aria-hidden="true">→</span></a>
</div>
</div>
<fieldset class="model-browser">
<legend class="visually-hidden">Preview a bundled model</legend>
<input class="model-choice" type="radio" name="preview-model" id="preview-chess" checked>
<label for="preview-chess">Chess</label>
<input class="model-choice" type="radio" name="preview-model" id="preview-head">
<label for="preview-head">Head</label>
<figure class="model-preview chess">
<img src="media/chess-viewport.webp" width="1664" height="1040" alt="Chessboard and ranks of light and dark NURBS pieces" fetchpriority="high">
<figcaption>Multipatch volumes</figcaption>
</figure>
<figure class="model-preview head">
<img src="media/head-viewport.webp" width="1664" height="1040" alt="Two heads modeled with surface patches, one wearing a cap and the other glasses">
<figcaption>Freeform surfaces</figcaption>
</figure>
</fieldset>
</div>

<h2 id="browse-the-documentation">Choose your workspace</h2>
<nav class="doc-paths" aria-label="Choose a documentation section">
<a href="page/fortran/index.html"><i class="fas fa-code" aria-hidden="true"></i><strong>Fortran API</strong><span>Build geometry in an fpm project, evaluate basis functions and derivatives, and export your results.</span></a>
<a href="page/how-to/index.html"><i class="fas fa-cube" aria-hidden="true"></i><strong>ForCAD Studio</strong><span>Edit geometry in the viewport, inspect connected patches, and generate meshes.</span></a>
</nav>

<h2 id="worked-problems">Worked problems</h2>
<ul>
<li><a href="page/fortran/project-a-point.html">Project a point onto a surface</a> to find its parameters and distance, or <a href="page/fortran/fit-a-curve.html">fit a curve to measured points</a>.</li>
<li><a href="page/fortran/curve-length-sensitivity.html">Differentiate curve length</a> with respect to a control height and its log weight.</li>
<li><a href="page/iga/solve-two-volumes-poisson.html">Solve Poisson’s equation across two volumes</a> and check the field and flux on their shared face.</li>
</ul>

<h2 id="library-and-examples">Reference and further reading</h2>
<div class="doc-shortcuts">
<a href="page/iga/index.html">Isogeometric analysis (IGA)</a>
<a href="page/reference/studio-fortran.html">Studio and Fortran equivalents</a>
<a href="page/reference/api.html">Public API reference</a>
<a href="page/reference/examples.html">Bundled examples</a>
<a href="page/explanation/index.html">Geometry concepts</a>
<a href="page/assistant/index.html">Studio AI Assistant</a>
</div>
