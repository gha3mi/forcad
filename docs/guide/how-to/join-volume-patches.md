title: Join two volume patches
category: Multipatch assembly

Create one volume assembly from two patches that share a face. To reproduce the example, [extrude a block](extrude-a-profile.html) and [split it](split-a-patch.html) halfway along U. Its lower and upper pieces share a face with aligned V and W directions.

1. Choose **Multipatch → New multipatch**. Set **Patch type** to **Volume**, name the assembly and check both volumes under **Members**.
2. Under **Interfaces**, choose **Add**. Select the lower piece as **Patch A** and its **U max** boundary. Select the upper piece as **Patch B** and its **U min** boundary.
3. For this split, keep **Parametric** continuity at **Order 0** and leave the orientation boxes clear. Choose **Check this interface**, then **Add to draft**.
4. Choose **Check interfaces**. When the row reads **Passed**, choose **Create multipatch**. Select the new entry under **Scene → Multipatches** to see the member and interface counts.

<figure class="doc-figure" markdown="1">

[![New multipatch dialog with two volume members selected and one C0 interface marked Passed between U max and U min](../../media/volume-interface.webp)](../../media/volume-interface.webp)

<figcaption>The selected row identifies the two faces. The detail below the table confirms their orientation and continuity.</figcaption>
</figure>

If the check fails, select the interface and choose **Edit**. Verify the two boundary labels and the orientation of their tangential directions before checking again. [How interface orientation works](../explanation/multipatch.html)

## Fortran equivalent

[Register and verify the two volume patches in Fortran](../fortran/join-volume-patches.html). To generate a conforming mesh in Studio, [mesh the assembly](mesh-a-multipatch.html).
