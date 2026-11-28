title: Connect an agent

1. Open a project. Use **Project → Folder** to check or set its working folder.
2. Open **Assistant** in the title bar. In **Agent**, choose an installed provider. ForCAD Studio can discover Codex CLI, OpenCode CLI, Claude Code CLI, and models served locally by Ollama.
3. Choose a **Model**. Where the provider exposes **Effort**, choose the reasoning effort for your request. If an agent or model is missing, check its installation, sign-in, or local service, then click the circular **Refresh available agents and models** button.
4. Click **Tools: Off** to switch to **Tools: Live**. Choose **Ask for approval** if you want to review project changes.
5. Ask “What geometry is in my open project?” and press **Enter**. An inspection in the activity stream confirms that the agent read the design. **Shift+Enter** inserts a line break.

<figure class="doc-figure" markdown="1">

[![Assistant controls with Codex CLI, GPT-6.2-Sol, live tools, and Ask for approval selected](../../media/assistant-agent.webp)](../../media/assistant-agent.webp)

<figcaption>Codex CLI with an OpenAI model that supports image input. The available models depend on the selected agent.</figcaption>
</figure>

A chat-only model can reply to text. To inspect or edit the open design, select a model that supports tool calling.
