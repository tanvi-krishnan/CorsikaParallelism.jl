# CORSIKA Change Log (mods vs baseline)

This file tracks modifications in:
- Modified source: `/n/holylfs05/LABS/arguelles_delgado_lab/Everyone/tkrishnan/corsika-mods`
- Baseline source: `/n/holylfs05/LABS/arguelles_delgado_lab/Everyone/tkrishnan/corsika-baseline-source`

## Current delta summary

### 1) `applications/c8_air_shower.cpp`
- Added `FirstInteractionStopper` process (file-scope class) to stop after first interaction.
- Added CLI option `--generation` (default `0`) for recursive cascade depth control.
- Added CLI flag `--stop-after-first-interaction`.
- Inserted `firstIntStopper` into the process sequence near the end.
- Wrapped `EAS.run()` with `try/catch` for `"FirstInteractionStop"`.
- In stop-after-first-interaction catch path, explicitly call `output.endOfShower()` before continuing so shower outputs are finalized/flushed.

### 2) `corsika/detail/modules/writers/InteractionWriter.inl`
- Extended interaction parquet rows to include:
  - `energy` (GeV, total energy)
  - `x`, `y`, `z` (m, interaction location in observation-plane coordinates)
- Extended parquet schema (`addField`) with matching columns: `energy`, `x`, `y`, `z`.
- Secondary energy call uses `particle.getEnergy()` (not `getTotalEnergy()`).

## Notes
- The explicit `output.endOfShower()` call was added to fix incomplete/empty interaction outputs when stopping via exception before normal `Cascade::run()` shower finalization.
- This file should be updated each time we modify `corsika-mods` relative to baseline.
