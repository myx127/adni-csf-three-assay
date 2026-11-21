# Paired three-assay CSF protein associations in Alzheimer’s disease

Code supporting **Same protein targets, different association patterns: a paired three-assay study of CSF protein associations with p-tau and longitudinal cognitive change across the Alzheimer’s disease spectrum**.

Version 1.0.0. The canonical repository is [https://github.com/myx127/adni-csf-three-assay](https://github.com/myx127/adni-csf-three-assay). Fixed releases are listed on the [Releases page](https://github.com/myx127/adni-csf-three-assay/releases). A Zenodo DOI will be added after archival.

## What is included

- The original analysis scripts from twelve study stages, with stage 11 containing the final 34-protein cognitive analysis and stage 12 the exploratory adjustment analyses.
- Stand-alone numerical kernels taken from the final pathology, clinical and paired-adjustment engines. Their function bodies are unchanged; data-dependent configuration imports are replaced with explicit numerical imports.
- Code and aggregate plotting inputs for all five main figures and sixteen supplementary figures.
- The 92 aggregate supplementary-table datasets and their field descriptions.
- Exact package versions from the environment used to verify this release, input-file release names, source hashes, and software-to-result crosswalks.
- A synthetic-data example and numerical/aggregate-result checks.

Participant-level data, participant identifiers, visit-level data, derived individual scores, fitted participant-level design objects, raw proteomic matrices and authentication material are not distributed.

## Start without ADNI data

Python 3.11 is recommended. From this directory:

```bash
python -m pip install -r requirements.txt
python examples/run_synthetic_demo.py
python -m pytest -q
python scripts/check_published_results.py
python scripts/render_figures.py --figures main1 main2 main3 main4 main5
```

The example is simulated and does not reproduce the article’s estimates. Figure rendering uses the already reported aggregate estimates; it does not refit models or recompute P values. English is the default. Chinese rendering requires a suitable CJK font.

## Analyses requiring authorized data

Read `docs/DATA_ACCESS.md`, `docs/REPRODUCIBILITY.md` and `metadata/workflow_catalog.csv`. The frozen workflows preserve the original stage names and numerical programs for traceability. They depend on governed source schemas, earlier-stage products and private checksums. Their verification-only wrappers are not a one-command pipeline from a fresh data download.

Use `scripts/prepare_private_workspace.py` to prepare a separate authorized workspace and local input manifests. Never run the frozen workflows in the author’s original analysis directories. Study-design calibration can require millions of draws; it is separate from the quick-start checks.

The current preparation verifies numerical kernels, synthetic execution, aggregate results and figure rendering. It does not rerun the entire participant-level analysis or all large calibration grids. Such reruns require the governed study releases and all earlier-stage derived objects described in the workflows.

## Interpretation

The final analyses distinguish shared association direction from equal association magnitude. Primary memory inference retains 68 hypotheses, including failed-calibration hypotheses represented by P=1. A non-significant joint memory test does not demonstrate absence of shared biological effects. NPTX2, TAFA5 and TIMP3 adjustment analyses are exploratory. NPTX1/KLK6 peptide evidence addresses a separate pathological-association branch.

## Citation and license

Software citation metadata are provided in `CITATION.cff` and `.zenodo.json`. The article has no assigned DOI in this package. The code is provided under the MIT License; ADNI data remain subject to ADNI access requirements.
