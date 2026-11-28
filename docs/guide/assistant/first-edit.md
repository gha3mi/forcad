title: Approve your first edit

Rotate the bundled Head model by 15 degrees, then restore its original orientation with Undo. You will need a [connected agent](connect-agent.html) with tool calling available.

1. On Home, open **Head** under **Samples**. Use **Project → Save as** to save a practice copy in your own project folder.
2. Open **Assistant** in the title bar. Enable **Tools: Live** and choose **Ask for approval**.
3. Send the following message.

   > Rotate the whole model 15 degrees about the world Z axis. Keep the camera in place.

4. Read the **ForCAD Studio approval** request. Check that the transform targets the whole scene, the axis is Z, and the angle is 15 degrees. Choose **Allow once**. If the proposed values differ, choose **Deny** and clarify your request.
5. Check the direction the two faces point relative to the viewport axes. Choose **Undo** in the title bar to restore the original orientation.

<figure class="doc-figure" markdown="1">

[![ForCAD Studio asks for approval to rotate the scene 15 degrees about Z, with Allow once and Deny controls](../../media/assistant-approval.webp)](../../media/assistant-approval.webp)

<figcaption>The approval prompt identifies the requested transform before ForCAD Studio changes the model.</figcaption>
</figure>
