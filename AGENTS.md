# AGENTS.md

This file is the standing instruction set for future coding agents working in this repository.

## Project Role

SEGD_reader is a reader/parser library, not a timing-correction engine. Its default job is to expose SEG-D header and trace contents faithfully, convert them into a practical Python/ObsPy shape, and avoid silent interpretation beyond the behavior already implemented in the package.

Do not add automatic correction, inversion, resampling, or timing repair logic unless the user explicitly requests a separate feature and the change is documented and tested as such.

## Current Scope

The package currently supports:

- SEG-D rev 3.0
- SEG-D rev 2.1

It has been tested on:

- Sercel: rev 3.0 and rev 2.1
- Stryde: rev 3.0
- SmartSolo: rev 2.1

Treat other manufacturers or revisions as potentially unsupported, even if the generic SEG-D structures appear similar.

## Repository Structure

- `SEGD_reader/segd_reader.py`
  Core parsing logic and `SEG_D_to_stream()`.
- `SEGD_reader/utils.py`
  Utility transforms, currently including Galperin-to-ENZ rotation.
- `SEGD_reader/json/rev_*/header_specs_rev_*.json`
  Versioned JSON header specs derived from the SEG-D documentation.
- `examples/data/`
  Bundled example files that should be preferred for lightweight regression tests.
- `segd_docs/`
  Source SEG-D PDFs. Treat these as the reference source when code or schema comments disagree.
- `tests/`
  Automated regression tests. Extend this rather than relying only on notebooks.

## Core Behavior and Assumptions

### Reader shape

- `SEG_D_Reader.get_general_headers_only()` reads the general headers only.
- `SEG_D_Reader.get_data_plus_key_headers()` returns a nested dictionary with:
  - file-level metadata under `file_header`
  - channel-set descriptions
  - per-trace header metadata
  - trace samples
- `SEG_D_to_stream()` converts parsed output into an ObsPy `Stream`.

### Version handling

- The reader discovers SEG-D version from General Header #2 unless a forced version is supplied.
- If an exact header-spec JSON is not available:
  - versions below 3 fall back to rev 2.1
  - versions 3 and above fall back to rev 3.0
- That fallback is permissive, not authoritative. When touching version logic, prefer preserving explicit warnings over silently broadening support claims.

### Time handling

- Rev 3.0 absolute record time is derived from General Header #3 `timeZero`, interpreted as SEG-D/GPS-epoch microseconds and converted to UTC using the built-in leap-second table.
- Rev 2.1 record time is reconstructed from year/day/hour/minute/second fields in General Header #1.
- `Trace.stats.starttime` currently comes from `file_header['record_timezero_utc']`.
- Optional per-trace Timestamp Header metadata may be exposed in parsed headers, but the default stream conversion does not currently recompute `Trace.stats.starttime` from it.

### Time Drift Header policy

- Expose raw timing metadata without silently altering record timing.
- Do not change `record_timezero_utc`, `Trace.stats.starttime`, sample spacing, or waveform alignment from Time Drift Header metadata alone.
- Time Drift Header values are provenance for downstream interpretation and comparison against independently measured drift.
- If future work adds optional timing correction, it must be explicit, documented, and tested separately from the default reader path.

For rev 3.0 Time Drift Header metadata, preserve the raw fields:

- `timeOfDeployment`
- `timeOfRetrieval`
- `timeOffsetDeployment`
- `timeOffsetRetrieval`
- `timedriftCorrected`
- `correctionMethod`
- `headerBlockType`

These should be exposed as:

- low-level parsed output: `trace_header['timeDriftBlock']`
- ObsPy stream output: `tr.stats.segd['timeDriftBlock']`

This metadata-exposure path has been tested on Stryde rev 3.0 data and should be treated as validated there before claiming support elsewhere.

### Calibration and scaling behavior

- `SEG_D_to_stream()` defaults to `use_descale_multiplier=True`.
- Descale multiplier application is already part of the package’s default behavior and changes trace amplitudes.
- Sensitivity correction is **not** applied by default; it only occurs when `apply_sensitivity_correction=True`.
- The stream metadata records whether descale/sensitivity corrections were applied via:
  - `tr.stats._descale_multiplier`
  - `tr.stats._sensitivity`
  - `tr.stats._units`
- Do not conflate timing metadata exposure with amplitude correction behavior when editing tests or docs.

### Station/component naming assumptions

- If `serial_to_station_name_dict` is provided, station/network/component naming is derived from that mapping.
- Without that mapping:
  - `station` defaults to the parsed `serialNumber` if available
  - otherwise it falls back to `line_point`
- Channel band codes are inferred from sample rate thresholds.
- Channel orientation suffixes come from `SENSOR_CODE` unless remapped through `serial_to_station_name_dict`.

### Merge and gap behavior

- `SEG_D_to_stream()` may concatenate multiple rev 2.1 trace segments for the same sensor before stream creation.
- `st.merge()` is attempted near the end of `SEG_D_to_stream()`.
- Merge errors are currently printed and ignored rather than raised.
- `remove_gaps=True` has two separate behaviors:
  - skips traces with invalid pre-1980 per-trace time-zero values
  - may fill missing 3C components with zero traces when `serial_to_station_name_dict` is provided
- `remove_stations_with_zero_data=True` removes stations whose traces are all zero or NaN before the optional component fill step.

### Integer conversion behavior

- `SEG_D_to_stream()` only converts samples to `int32` when all trace values are integral after all enabled corrections.
- Tests that compare waveform data should account for the possibility of float output depending on descale/sensitivity paths.

## Manufacturer-Specific Assumptions

- Rev 3.0 is mostly driven by generic JSON header specs.
- Rev 2.1 still contains manufacturer-specific assumptions in code.
- At present, additional rev 2.1 trace-extension parsing for latitude/longitude/elevation/serial number/sensitivity is only implemented for Sercel (`manufacturersCode == 13`).
- README already notes that sensitivity and location info for rev 2.1 are manufacturer-defined and only known here for Sercel. Preserve that warning unless the user has verified another layout.

Do not generalize rev 2.1 manufacturer-specific parsing from one vendor to all vendors without evidence and tests.

## Rev 3.0 Header Mapping

- `0x42`: Timestamp Header
- `0x43`: Sensor Calibration Header
- `0x44`: Time Drift Header

When schema text, code comments, and `segd_docs` disagree, treat the PDFs in `segd_docs/` as the source of truth and update the derived schema/comments to match.

## Testing Expectations

- Use test-first changes for parser behavior.
- Add regression coverage for bundled example files whenever practical.
- Prefer `examples/data/1.101_19_00_00_SN_2402041994_0.rsamp.segd` for rev 3.0 time-drift regression coverage unless a better bundled fixture is added.
- Verify metadata exposure separately from waveform-transform behavior.
- Guard tests should prove that metadata-only changes do not silently change timing or trace samples.
- When changing docs or schema wording around block semantics, add or update a regression check if the wording is safety-relevant.

## Documentation Expectations

- Keep comments, README notes, JSON schema text, and code behavior aligned with the SEG-D documents in `segd_docs/`.
- Be explicit about whether behavior is spec-derived or manufacturer-specific.
- If a field is exposed but not applied, say that plainly in the docs.
- If a behavior is approximate, permissive, or fallback-based, document that instead of implying strong correctness.

## Editing Guidance

- Prefer small, local changes in `segd_reader.py`; this file already mixes format parsing, metadata shaping, and stream conversion logic.
- Avoid broad refactors unless the user explicitly asks for cleanup.
- Preserve existing public function names and default argument behavior unless the user requests an API change.
- If you change any semantics in `SEG_D_to_stream()`, review README and tests in the same change.

## Verification Checklist

Before wrapping up a parser change, verify at minimum:

- the targeted pytest regression passes
- metadata lands in both the low-level parsed structure and stream metadata when intended
- timing fields remain unchanged when the change is meant to be metadata-only
- README and schema wording do not contradict `segd_docs`
