# Changelog

## version 4.1.0 - 10/02/2026
### Summary
New porosity model (DDP)

### Details
- add depth dependent porosity model (DDP)
- improve flux computation with phi_min value
- more detailed documented tests
- bug fixes

## version 4.0.2 - 10/11/2025
### Summary
Unit tests and complete doc

### Details
- add unit tests (CI runs) for 20% of the common routines
- add merge and validation tests for models
- CI builds all docs (user and developer), including HTML
- bug fixes

## version 4.0.1 - 30/10/2025
### Summary
Small fixes (see details in CHANGES.md)

### Details
- add test pipeline on public repository
- add unit tests for geometry routines
- add changes.md
- doc correction
- python scripts revisited

## version 4.0.0 - 17/10/2025

### Summary
New open source release including porosity models

### Details
- sw2dSolver: core solver of the software, including
	- SWES: traditional Saint-Venant equations
	- LISWES: local inertia approximation to the shallow water equations
	- SP: single porosity shallow water equations
	- DIP: dual integral porosity shallow water equations
	- available source terms: rain, infiltration, wind, friction
	- available schemes: Godunov and MUSCL-EVR
- sw2dConverter: converting 2dm mesh into binary information for the solver
- sw2dPorosityMapper: tool to build porosity map for the SP and DIP models
- sw2dPostProcess: post-processing tools (computing max values, wetting times, flood duration, etc.)
- sw2dZone: tool to build heterogeneous fiels for source terms (friction, wind, etc.)
- documentation (pdf and html) for users and developers
- test cases library