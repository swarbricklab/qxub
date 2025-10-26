# qxub Documentation Guide

This guide defines how we organize, write, and validate documentation across the project.
It also establishes conventions for runnable examples and the “Documentation Chat Mode”.

## Principles

* Be lean and user-focused (80/20 rule).
* Avoid duplication; prefer a single source of truth and link to it.
* Keep quick start in the top-level `README.md` under ~60 lines.
* Reflect reality: docs must match the current implementation and CLI help.
* Runnable examples must be safe by default and easy to test automatically.

## Documentation Structure

* Top-level `README.md` (quick start)
  * What qxub is, how to install, the most common one-liner, and links to deeper docs.

* `docs/` (user docs)
  * `README.md` – index of user docs.
  * `examples.md` – the authoritative collection of common, runnable examples.
  * `configuration.md` – config precedence, XDG paths, and essential options only.
  * `aliases.md` – alias usage patterns.
  * `platform_configuration.md` – platform definitions and queue selection rules.
  * `remote-execution.md` – SSH/remote execution usage (emerging feature).
  * `option-placement.md` – unified CLI structure and option placement.
  * `shortcuts.md`, `smart-quotes.md` – focused user topics.
  * `tutorial/` – step-by-step guides (can include runnable snippets when feasible).

* `docs/dev/` (developer docs)
  * Architecture and design references (schemas, threading, package structure).
  * Non-runnable deep dives kept concise; update when implementation changes.

* `docs/RELEASE_NOTES_*.md` (release notes)

* `examples/` (worked scenarios)
  * Optional scripts that demonstrate end-to-end flows; can be referenced by docs.

* Docstrings
  * Keep short, precise, and include doctest examples when helpful.

## Runnable Examples Convention

We automatically extract and test runnable examples from Markdown files.

* Mark runnable code blocks with a fenced header containing both the shell language and the `runnable` tag:

  ```
  ```bash runnable
  qxub --help
  ```
  ```

* Allowed languages: `bash`, `sh`.
* Optional tags: Add a first-line comment to categorize blocks, e.g.:
  * `# tags: hpc` – requires access to an HPC environment; skipped by default.
  * `# tags: slow` – long-running; skipped by default.
  * `# tags: destructive` – contains potentially destructive operations; rejected.

Safety rules enforced by the test runner:
* Blocks containing obviously destructive commands are rejected (e.g., `rm -rf`, `sudo`, `mkfs`, `dd if=`).
* By default, blocks tagged `hpc` or `slow` are listed but skipped unless explicitly enabled.
* Prefer dry-run flags in examples where available.

Scope for extraction:
* Top-level `README.md` and files under `docs/` (including `docs/tutorial/`).
* Developer docs under `docs/dev/` are not scanned by default.

## Validation Tooling

Run the example tests with:

```bash
tests/test_docs_examples.sh --list       # enumerate runnable snippets
tests/test_docs_examples.sh --run        # run safe, untagged snippets
tests/test_docs_examples.sh --run --include-tag hpc  # opt-in categories
```

Guidelines for authors:
* Add runnable examples to `docs/examples.md` when feasible.
* Keep examples self-contained and short.
* Prefer safe defaults and dry runs.
* After updating docs or CLI behavior, run the example tests.

## Update Checklist (for PRs)

When adding or changing functionality:
1. Update user docs (examples, options, configuration) as needed.
2. Update tutorials only if user flows change.
3. Update dev docs if architecture or contracts change.
4. Add/adjust runnable examples and run `tests/test_docs_examples.sh`.
5. If docstrings changed, ensure doctests still pass (when enabled in CI).

## Documentation Chat Mode Contract

When using Documentation Chat Mode:
* Always consult this guide first to determine where content belongs.
* Keep responses lean and consistent with the structure above.
* Validate claims against:
  * CLI help (`qxub --help`, subcommand `--help` where applicable), and
  * The files under `qxub/` as the source of truth.
* Prefer updating `docs/examples.md` for examples.
* For runnable examples, include the `bash runnable` fence and appropriate tags.
* Suggest running `tests/test_docs_examples.sh` to verify examples.

Disallowed patterns:
* Duplicating content across multiple pages without a clear reason.
* Including long implementation details in user docs (link to dev docs instead).

This guide is the canonical reference for documentation organization, conventions, and validation in qxub.
