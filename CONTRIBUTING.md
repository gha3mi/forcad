# Contributing to ForCAD

Contributions to ForCAD are welcome and thanks for taking the time to help improve the project!

## Table of Contents
- [Contributing to ForCAD](#contributing-to-forcad)
  - [Table of Contents](#table-of-contents)
  - [How to Contribute?](#how-to-contribute)
    - [Reporting Issues](#reporting-issues)
    - [Suggesting Enhancements](#suggesting-enhancements)
    - [Implementing Changes and Adding Features](#implementing-changes-and-adding-features)
    - [Sharing Examples](#sharing-examples)
    - [Release version](#release-version)
  - [Styleguides](#styleguides)
    - [Code Style](#code-style)
    - [Documentation](#documentation)
    - [Git Commit Messages](#git-commit-messages)

## How to Contribute?

### Reporting Issues

If you encounter any problems, please [open an issue](https://github.com/gha3mi/forcad/issues) on GitHub. Be sure to include detailed information about the issue.

### Suggesting Enhancements

Do you have ideas for improving ForCAD? Please [open an issue](https://github.com/gha3mi/forcad/issues) on GitHub and describe your proposed enhancements in detail.

### Implementing Changes and Adding Features

If you'd like to add new features, fix bugs, or enhance existing functionality, consider submitting a pull request (PR). Here's how to proceed:

1. **Fork** the repository and create a new branch from `main` to work on your changes.
2. Implement your modifications and enhancements.
3. Ensure that your code adheres to the project's coding style and guidelines.
4. Write clear and descriptive commit messages.
5. Submit your PR with an explanation of the changes introduced.

### Sharing Examples

Contribute to ForCAD by submitting a pull request (PR) with examples. Ensure that your examples are documented.

### Release version

The root [VERSION](VERSION) file supplies the package versions, the Studio About dialog, Assistant handshakes, and the documentation footer. To select a release number, edit that file once and rebuild. `bash release.sh --dry-run` previews the release notes. When the file still matches the last release tag, the helper proposes a version based on the intervening commits.

## Styleguides

### Code Style

- Maintain consistency with the existing code style and structure.
- Write clear, concise, and well-commented code.
- Thoroughly test your changes before submission.

### Documentation

- Choose the Diátaxis section under `docs/guide/` that matches the reader's goal. Use tutorials for guided learning, how-to guides for tasks, reference for lookup, and explanation for concepts. [Diátaxis](https://diataxis.fr/)
- Use sentence-case headings, exact Studio labels, runnable code, and expected results. Give screenshots meaningful alt text and keep them in `docs/images/`. [Google developer style guide](https://developers.google.com/style/)
- Document interfaces exported by `use forcad` in the Fortran source. When a helper needs public visibility between core modules only, keep its access declaration in the existing `#ifndef FORCAD_PUBLIC_DOCS` group. Use LaTeX for mathematical expressions and code formatting for Fortran identifiers. Where applicable, include array shapes, parameter intervals, units, and error behavior. [FORD documentation syntax](https://forddocs.readthedocs.io/en/stable/user_guide/writing_documentation.html)
- Give each how-to guide one task and an outcome the reader can check. Put it under ForCAD Studio (`how-to/`), `fortran/`, or `assistant/`. Keep runnable public API examples in `fortran/` and link the equivalent Studio steps in both directions. In the comparison reference, use a dash for a Studio action without a direct public API.
- Capture the relevant control in a clean practice project. Exclude accounts, personal paths, and private chat content. When a model selector is relevant, show an OpenAI model. Keep the visible control labels in the instructions.
- Add each new guide to its section's `index.md` and link it from a related task. Keep published page paths and heading anchors usable when reorganizing topics.
- Run the steps or compile the examples you change. Build with `pixi run docs` and open `docs/ford/index.html`. Edit the landing page in `docs/index.md`, layouts in `docs/templates/`, and colors and spacing in `docs/ford.css`.
- For layout or theme changes, check a guide, a public API page, and search. Try the keyboard, Light and Dark, and a 320-pixel viewport. [WCAG reflow](https://www.w3.org/WAI/WCAG22/Understanding/reflow.html)

### Git Commit Messages

Use the [Conventional Commits](https://www.conventionalcommits.org/) specification to format your commit messages. This helps keep the project history readable and enables automatic changelog generation and versioning.

**Format:**

```
<type>(optional-scope): <description>

[optional body]
```
**Guidelines:**

* Keep the subject line **under 72 characters**.
* Use lowercase for `<type>` and `<scope>`.

**Common `<type>` values:**

  * **feat**: a new feature
  * **fix**: a bug fix
  * **refactor**: code change that neither fixes a bug nor adds a feature
  * **docs**: documentation only changes
  * **style**: changes that do not affect the meaning of the code (white-space, formatting, etc.)
  * **perf**: a code change that improves performance
  * **chore**: other changes that don't modify `src`, `example` or `test` files
