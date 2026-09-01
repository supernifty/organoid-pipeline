# Repository instructions

## Spec-driven development

- Before implementing any feature, fix, refactor, or other code change, create or update a written specification in `specs/`.
- Treat the specification as the source of truth for the work. It must state the problem, goals, scope, acceptance criteria, and verification approach clearly enough to guide implementation.
- Do not begin implementation until the relevant specification exists.
- Keep the specification current when requirements, scope, or implementation decisions change.

## Temporary files

- Create task-specific temporary files and directories inside the current workspace, preferably under `tmp/codex/`.
- Never create task artifacts under `/tmp`, `/private/tmp`, or the system temporary directory unless explicitly approved.
- Clean up temporary artifacts after verification.
- Before copying large files or directories, estimate their size and available workspace capacity.

## Mathematics in CLI responses

- Do not use LaTeX delimiters or LaTeX commands in CLI-facing responses.
- Write mathematics as readable plain text, using Unicode symbols such as −, ×, ÷, ≤, ≥, ≈, →, and Σ where they improve clarity.
- Put multi-line equations in fenced plain-text code blocks so spacing is preserved.
- Define symbols in prose and prefer forms such as `d = N_ID / (N_ID + N_SBS)`.
- Use an ASCII fallback when a Unicode symbol may be ambiguous or unsupported.
