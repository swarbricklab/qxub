---
description: 'Chat mode for maintaining documentation quality'
tools:
  # Core, read-only
  - search
  - codebase
  - fileSearch
  - readFile
  - editFiles
  - usages
  - problems
  - changes
  - fetch
model: GPT-5
---
# Documentation Chat Mode

This chat mode focuses on maintaining clear, consistent, and validated documentation.

## Operating Rules

1) Always load and follow `docs/documentation-guide.md`.
2) Keep user docs lean (80/20) and avoid duplication; link to a single source of truth.
3) Validate statements against the implementation:
   - Prefer reading CLI help (e.g., `qxub --help` and subcommand `--help`).
   - Consult code under `qxub/` for authoritative behavior.
4) Prefer adding runnable examples to `docs/examples.md` using the `bash runnable` fence.
5) After editing examples, ask the user to run:
   - `tests/test_docs_examples.sh --list` then `--run` to validate snippets.
6) For major behavior changes, remind to update release notes and relevant tutorials.
7) Avoid documentation bloat:
   - Update an existing document where possible rather than creating a new file.
   - Only create new docs when explicitly requested or when a clear gap exists.
   - If no obvious home, propose a target and seek clarification before proceeding.

## Scope and Placement

- Quick starts stay in the top-level `README.md` (under ~60 lines) and link into `docs/`.
- User guides live in `docs/*.md` and `docs/tutorial/`. Dev architecture stays in `docs/dev/`.
- Keep implementation detail out of user docs; summarize and link to dev docs when needed.

## Runnable Examples Convention

- Mark runnable shell examples as:

  ```
  ```bash runnable
  # tags: hpc
  qxub --help
  ```
  ```

- Use tags (`hpc`, `slow`, etc.) to control the test runner.
- Ensure commands are safe (prefer dry-run flags) and self-contained.

## Consistency Checklist

When updating docs for a feature:
- [ ] Does `qxub --help` reflect the documented options and subcommands?
- [ ] Are examples updated in `docs/examples.md` and correctly marked as runnable?
- [ ] Are tutorials accurate for the updated flows?
- [ ] Are dev docs updated if contracts/architecture changed?
- [ ] Have you run `tests/test_docs_examples.sh` to validate runnable examples?

## Authoring Notes

- Use short paragraphs and code examples over long prose.
- Prefer concrete commands that users can copy and try.
- If the same topic appears in multiple places, consolidate to one page and link.
