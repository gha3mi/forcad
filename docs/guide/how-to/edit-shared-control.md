title: Edit a shared multipatch control
category: Patch editing

Use the assembly's control editor to change an interface point or weight across its connected patches.

1. Select the assembly row under **Scene → Multipatches**. Open **Multipatch → Points + weights**.
2. Locate the **Control ID** and double-click the coordinate or **Weight** you want to edit. Enter the new value.
3. If the assembly has a \(G^1\) join, choose **Preview**. The table identifies each dependent control, the affected field, and its **Before** and **After** values. Preview again if you edit a value.
4. Choose **Apply** and inspect the joined patches. **Undo** restores the edit and the dependent changes together.

See [Multipatch interfaces](../explanation/multipatch.html) for how the interface relates controls across the join.
