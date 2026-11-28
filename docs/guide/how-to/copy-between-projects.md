title: Copy objects between project tabs
category: Projects and exchange

1. Open both `.fcad` files with **Project → Open**. Each project has its own tab.
2. In the source tab, select the geometry in **Scene**. Use **Ctrl+click** to add or remove objects from the selection. Choose **Project → Copy** (**Ctrl+C**).
3. Click the destination project’s tab and choose **Project → Paste** (**Ctrl+V**). The pasted objects are selected in the destination scene.
4. Save the destination project (**Ctrl+S**).

<figure class="doc-figure">
  <a href="../../media/paste-project-result.webp">
    <picture>
      <source media="(max-width: 32rem)" srcset="../../media/paste-project-tabs.webp" width="307" height="33">
      <img src="../../media/paste-project-result.webp" width="478" height="348" alt="After copying a plate from Profiles.fcad, Layout.fcad is the active tab and has an asterisk.">
    </picture>
  </a>
  <figcaption>Plate is selected under Surfaces in Layout.fcad. The asterisk marks the destination’s unsaved change.</figcaption>
</figure>

To transfer an entire model, [import the FCAD project](import-fcad.html).
