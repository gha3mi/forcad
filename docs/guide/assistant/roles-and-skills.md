title: Agent roles and skills

The **Agent** selector chooses the provider that answers your request. A **role** gives that agent a responsibility. A **skill** is a workflow the agent can consult for a particular modeling task. You can use skills with a single agent or a team. In the checkout, role briefs live in `assistant/agents/` and each skill has an `assistant/skills/<id>/SKILL.md` file.

## Who does what?

With **Tools: Live**, the Assistant receives the shared [Studio instructions](https://github.com/gha3mi/forcad/blob/dev/assistant/AGENTS.md). They direct it to check the active tab, units, selection, and object IDs, preview edits where available, follow the chosen approval policy, and inspect the result.

| Role brief | Responsibility |
| --- | --- |
| [Lead](https://github.com/gha3mi/forcad/blob/dev/assistant/agents/lead.md) | Checks the live design, considers the team’s findings, carries out edits in order, and reports what changed and what you can undo. |
| [Geometry reviewer](https://github.com/gha3mi/forcad/blob/dev/assistant/agents/geometry-reviewer.md) | Reviews geometry, dimensions, and interfaces. Flags missing evidence or continuity risks and proposes checks. |
| [Modeling planner](https://github.com/gha3mi/forcad/blob/dev/assistant/agents/modeling-planner.md) | Proposes construction steps, prerequisites, patch layout, and verification. |

In [Team mode](configure-team.html), the reviewer and planner work in parallel and return proposals to the lead. The lead makes the edits under your approval settings.

## Which skill fits my task?

| Skill | When to use it |
| --- | --- |
| [`scene-inspection`](https://github.com/gha3mi/forcad/blob/dev/assistant/skills/scene-inspection/SKILL.md) | Identify objects or mesh cells, measure geometry, or check the current selection and camera. |
| [`geometry-construction`](https://github.com/gha3mi/forcad/blob/dev/assistant/skills/geometry-construction/SKILL.md) | Build curves, surfaces, or volumes from dimensions or a reference. Identify any inferred dimensions. |
| [`patch-editing`](https://github.com/gha3mi/forcad/blob/dev/assistant/skills/patch-editing/SKILL.md) | Edit controls, weights, knots, degrees, or parameter ranges and preview dependent changes. |
| [`multipatch`](https://github.com/gha3mi/forcad/blob/dev/assistant/skills/multipatch/SKILL.md) | Check face orientation, join members, or refine a connected assembly. |
| [`structured-meshing`](https://github.com/gha3mi/forcad/blob/dev/assistant/skills/structured-meshing/SKILL.md) | Preview cell counts, mesh by target size, or check that a hex mesh conforms across interfaces. |
| [`project-exchange`](https://github.com/gha3mi/forcad/blob/dev/assistant/skills/project-exchange/SKILL.md) | Import geometry, work across project tabs, save or export a model, or capture the viewport. |

For example, select a volume assembly, enable **Tools: Live**, and ask “Use the structured-meshing skill to preview cell counts for this assembly. Show me which interface directions must match before generating the mesh.”

To see which skills your copy of Studio provides, ask the Assistant to list its skills. The read-only `studio_skill` tool lists available guides when called with `{}`. With `{"id":"structured-meshing"}` it returns that guide’s text. Reading a skill leaves the model untouched.
