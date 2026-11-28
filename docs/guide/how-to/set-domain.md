title: Set a parameter domain
category: Patch editing

Change the parameter interval of one direction while preserving the physical shape of the patch.

1. Select one patch and open **Patch → Set domain**.
2. Choose a direction and read **Current range**.
3. Enter **New lower** and **New upper**. Both bounds must be finite and the new upper bound must exceed the lower bound.
4. Choose **Set domain**. Reopen the dialog and select the same direction to verify the new range.

<figure class="doc-figure">
  <a href="../../media/set-domain.webp"><img src="../../media/set-domain.webp" width="450" height="240" alt="Set domain for a unit cube. U has a current range of 0 to 1 and the new bounds are 2 and 6."></a>
  <figcaption>Read the existing interval under <strong>Current range</strong> and enter its replacement in the two fields below.</figcaption>
</figure>

For example, changing the U interval from \([0,1]\) to \([2,6]\) changes the midpoint parameter from 0.5 to 4. Both values evaluate to the same physical point.

**Fortran equivalent.** [Change a parameter interval and evaluate the same physical point](../fortran/set-domain.html).
