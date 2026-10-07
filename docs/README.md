# Mathematical documentation

These documents provide a maintained home for mathematical explanations, definitions, derivations,
and useful exploratory reasoning in SSE, initially focused on QuaSSE. Coverage is selective: this is
not a complete specification or a certification of the implementation.

## Reading and placement

- [QuaSSE guide](quasse/guide.tex): accessible explanations of the model and its interpretation.
- [QuaSSE theory](quasse/theory.tex): precise conventions, assumptions, derivations, and proofs.
- [QuaSSE explorations](quasse/explorations.tex): developing ideas, alternatives, and unresolved issues.
- Existing [examples](../examples/QuaSSE_fossils_fixed_tree.md) retain commands and XML instructions;
  existing [validation reports](../validation/) retain detailed numerical studies and measurements.

Organize by subject, not file format. Markdown tutorials may live alongside the LaTeX files. Keep
small documents in one file; split topics into included files when their length makes this useful.
A future `literature/quasse-and-sse.tex` can discuss applications and relationships among published
methods once there is enough material. Specific sources supporting a derivation belong beside it,
using the shared `references.bib`, even if a literature survey also discusses that paper.

## Building

Install a TeX distribution with `pdflatex`, BibTeX, and `latexmk`. The documents use standard packages
including `amsmath`, `amssymb`, `geometry`, `microtype`, `booktabs`, and `hyperref`.

```sh
make -C docs
make -C docs theory       # also: guide, explorations
```

PDFs and auxiliary files go into the ignored `build/docs/` directory at the repository root.
Each document has its own auxiliary directory. Run `make -C docs clean` to remove only these
three documents' generated output. Keep sources and original figures in version control, not PDFs
or TeX auxiliary files. The PDFs' repository links assume they remain in `build/docs/`; distributing
PDFs separately may require providing the linked material as well.

## Writing and maintaining

Consult the relevant explanation when it helps development. When an authorized implementation
requires an existing explanation to change to remain accurate, make focused edits to that explanation.
This is a flexible resource: no documentation edit is required for every code change, and there is
no mandatory claim registry, implementation map, provenance ledger, or review sequence.

- Develop a coherent account of a topic rather than appending dated updates or changelog entries.
- Define notation and assumptions before using them. State time direction and parameter units.
- Distinguish identities and proofs from approximations, empirical checks, and tentative ideas.
  Agreement between implementations is evidence, not a proof of the mathematical model.
- Explain whether a method is implemented or proposed when that distinction affects understanding.
- Cite borrowed results and, where useful, their equation or section. Read the relevant source and
  explain differences in conventions; do not invent references or imply unperformed verification.
- Keep one detailed treatment of a result; explain it briefly elsewhere and point to that treatment.
- Use ignored locations for scratch material and report findings in the conversation.
  Do not track scratch material or move it into a tracked directory without explicit authorization.
- Include broader SSE mathematics only when it serves a concrete purpose for this collection.
- Use stable ordinary LaTeX labels for results worth referencing. Share notation and formatting only
  where actually useful; do not introduce a document-generation framework.

The guide favors interpretation, the theory document supplies justification, and explorations allow
unfinished reasoning. None should become an inventory of all possible future topics. Describe
limitations near the relevant claim without repeating disclaimers throughout the collection.

## Authorization for permanent records

Do not create or expand permanent development records unless the user explicitly requests them.
This includes validation reports, benchmark summaries, audit reports, session notes, decision
histories, and changelog entries. Authorization to investigate, implement, test, or benchmark does
not itself authorize preserving a report in the repository. Existing reports, directories, and
previous exceptions do not establish continuing permission.

Focused edits needed to keep existing user or mathematical explanations accurate are distinct from
appending a record of work performed. If a permanent record would be valuable, propose its purpose
and destination; do not create it pending approval. Placement guidance in this README describes
where authorized material belongs, not permission to create it.

## Exploration lifecycle

When the user explicitly requests creating or revising exploration documentation, an active
exploration should explain its question, motivation, current formulation, evidence,
alternatives, and unresolved issues. Use these as writing guidance, not a compulsory form with empty
fields. State its status in prose, such as experimental or adopted but not yet implemented.

When revising an authorized exploration after its question is settled, move established mathematics
into the theory document and user-facing
interpretation into the guide. Condense the exploration into the decision, important alternatives,
reasons for the choice, and remaining limitations. A brief history is useful when it explains how
reasoning changed. Preserve detailed rejected arguments only when they contain a useful result,
explain a subtle failure, or remain plausible options to revisit. Version control retains superseded
working detail. Do not casually discard existing evidence needed to understand a retained scientific
conclusion. This preservation guidance does not authorize creating additional reports or histories.

These conventions apply to this collection; existing validation reports need not be retroactively
rewritten to match its format.
