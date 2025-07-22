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

## Styleguides

### Code Style

- Maintain consistency with the existing code style and structure.
- Write clear, concise, and well-commented code.
- Thoroughly test your changes before submission.

### Documentation

- Follow the [FORD documentation styles](https://forddocs.readthedocs.io/en/latest/user_guide/writing_documentation.html) for documentation.

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