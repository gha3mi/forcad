title: Where Assistant data goes

The project working folder keeps the saved conversation and its references alongside your design. Choose **Project → Folder** to check which folder is active.

| Item | Where ForCAD Studio keeps it |
| --- | --- |
| Model and meshes | The project’s `.fcad` file |
| Explicitly saved chats | `.forcad/chats.ini` in the working folder |
| Copies of attached files | `.forcad/references/` in the working folder |
| Captured viewport images | `.forcad/screenshots/` in the working folder |

Sharing an entire working folder can include its hidden `.forcad` directory. To share just the saved model, send its `.fcad` file.

The selected agent receives the prompt and the context supplied for the request. Live tools can return information about the open design. A viewport capture excludes the Assistant conversation. Each provider controls how it processes or retains the material it receives. Deleting a locally saved chat removes its entry in ForCAD Studio.
