# Input schemas and scales

Official file releases are listed in `metadata/input_files.csv`; the original readers are retained in `workflows/`. NULISA protein records, SomaScan sample/analyte identifiers and TMT-MS sample/protein metadata are mapped to participant and CSF collection date. The earliest common date defines the paired index; matching does not establish an identical aliquot or assay-run date.

The clinical-engine frame includes `age_month_precision`, `sex_code2`, `education`, `apoe4_carrier`, `P_LOG`, `NULISA_PlateID`, `SOMASCAN_PlateId`, `TMT_plex` and `ptau_batch`. Longitudinal construction additionally requires actual cognitive dates, ADNI memory/executive composite scores and within-participant time information. The final templates include A (coefficient operators), K (stacked technical-cluster score operators), sigma (fixed protein scales), and Uv (scale-estimation influence scores).

The paired-adjustment template combines model-specific coefficients and scores and a contrast matrix L. Its covariance must include contributions shared between the models. Independent standard errors or differences between confidence-interval endpoints cannot replace that paired covariance.

The synthetic example supplies known fixed protein scales and sets Uv to zero for that example only; production source code retains estimated-scale uncertainty. It generates summarized annual memory changes, not a substitute for reconstructing actual ADNI longitudinal design templates.
