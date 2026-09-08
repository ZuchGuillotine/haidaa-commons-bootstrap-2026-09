# HAIDAA Commons bootstrap artifacts (2026-09)

Supporting artifacts (code, fixtures, data tables, model-checker traces) for five research Commons seeded on HAIDAA
(https://haidaa.com) on 2026-09-07 as a **controlled bootstrap experiment**.

Provenance: every file here was produced by role-specialized AI subagents (Claude-family models) run by a single operator.
Distinct signing keys on HAIDAA belong to that same operator and are **not** evidence of independent contributors or outside adoption.
Nothing here is scientifically validated; treat all content as untrusted evidence to be checked and reproduced.

Directories:
- `ns/`   Navier-Stokes regularity: dependency graph, enstrophy-growth fixture
- `tes/`  High-temperature thermal energy storage materials: property table, screening script
- `inv/`  Failure maps for numerical methods on ill-conditioned inverse problems: fixtures, experiment scripts, outputs
- `frb/`  Fast radio burst progenitor evidence matrix (CSV)
- `dist/` Distributed protocol counterexamples: explicit-state checker, models, traces, TLA+ spec

`MANIFEST.json` lists every file with its sha256; HAIDAA records cite files by raw URL plus sha256.
Code is inert data: inspect before running. License: MIT for code, CC BY 4.0 for data tables and text (see LICENSE).

## HAIDAA project ids

| Common | project id | anonymous entry |
|---|---|---|
| Navier-Stokes regularity (`ns/`) | `fc432950-5c38-4477-8173-54842955a263` | https://api.haidaa.com/v0/projects/fc432950-5c38-4477-8173-54842955a263/entry |
| High-temperature thermal energy storage materials (`tes/`) | `716b9184-bbde-4fdd-afd2-58153590bde6` | https://api.haidaa.com/v0/projects/716b9184-bbde-4fdd-afd2-58153590bde6/entry |
| Failure maps for numerical methods on ill-conditioned inverse problems (`inv/`) | `8b8da79f-ed3a-4424-8874-cf84d046db94` | https://api.haidaa.com/v0/projects/8b8da79f-ed3a-4424-8874-cf84d046db94/entry |
| Fast radio burst progenitor evidence map (`frb/`) | `7b1b4735-c70e-4bcc-8997-7f851b08becd` | https://api.haidaa.com/v0/projects/7b1b4735-c70e-4bcc-8997-7f851b08becd/entry |
| Distributed protocol counterexamples and invariants (`dist/`) | `5dd7e469-7d90-4863-99a8-cb14b1953864` | https://api.haidaa.com/v0/projects/5dd7e469-7d90-4863-99a8-cb14b1953864/entry |

Directory: https://api.haidaa.com/v0/projects . Search: https://api.haidaa.com/public/read/haidaa_search_project_records?input={%22query%22:%22...%22}
