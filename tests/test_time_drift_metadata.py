from pathlib import Path

import json
import numpy as np

from SEGD_reader import SEG_D_Reader, SEG_D_to_stream


REPO_ROOT = Path(__file__).resolve().parents[1]
REV3_SAMPLE = REPO_ROOT / "examples" / "data" / "1.101_19_00_00_SN_2402041994_0.rsamp.segd"
HEADER_SPECS = REPO_ROOT / "SEGD_reader" / "json" / "rev_3_0" / "header_specs_rev_3_0.json"
README = REPO_ROOT / "README.md"
AGENTS = REPO_ROOT / "AGENTS.md"

EXPECTED_TIME_DRIFT_KEYS = {
    "timeOfDeployment",
    "timeOfRetrieval",
    "timeOffsetDeployment",
    "timeOffsetRetrieval",
    "timedriftCorrected",
    "correctionMethod",
    "headerBlockType",
}


def _first_trace_entry(data):
    for channel_set, payload in data.items():
        if channel_set == "file_header":
            continue
        for line_data in payload.get("traceData", {}).values():
            for trace_data in line_data.values():
                return trace_data
    raise AssertionError("No trace data found in parsed SEG-D output")


def _read_raw_time_drift_block(path: Path):
    reader = SEG_D_Reader(str(path), verbose=False)
    reader.open_file()
    try:
        # Load header specs for this file version before using low-level helpers.
        reader.get_data_plus_key_headers()
        general_header3 = reader.read_general_header3()
        header_size = general_header3["headerSize"]
        start_byte = header_size

        channel_sets = []
        header_offset = 3 * 32
        while header_offset < header_size:
            header = reader.read_other_header(header_offset)
            if (
                isinstance(header, dict)
                and header.get("headerBlockType") == 48
                and header.get("headerBlockType2") == 49
                and header.get("headerBlockType3") == 50
            ):
                channel_sets.append(header)
            header_offset += 32

        bytes_per_sample = 3
        if reader.read_general_header1()["formatCode"] == "8058":
            bytes_per_sample = 4

        for description in channel_sets:
            for _ in range(description["numberOfChannelsThisSet"]):
                start_byte += 20
                for extension_index in range(description["numberTraceHeaderExtensions"]):
                    parsed_header = (
                        reader.read_trace_header_extension(start_byte)
                        if extension_index == 0
                        else reader.read_other_header(start_byte)
                    )
                    if (
                        isinstance(parsed_header, dict)
                        and parsed_header.get("headerBlockType") == 68
                    ):
                        return parsed_header
                    start_byte += 32
                start_byte += description["samplesPerChannel"] * bytes_per_sample
    finally:
        reader.close_file()

    raise AssertionError("Expected raw 0x44 Time Drift Header in rev-3 sample file")


def test_rev3_sample_contains_raw_time_drift_header():
    raw_block = _read_raw_time_drift_block(REV3_SAMPLE)

    assert set(raw_block) == EXPECTED_TIME_DRIFT_KEYS
    assert raw_block["headerBlockType"] == 68


def test_get_data_plus_key_headers_exposes_time_drift_block():
    reader = SEG_D_Reader(str(REV3_SAMPLE), verbose=False)
    reader.open_file()
    try:
        data = reader.get_data_plus_key_headers()
    finally:
        reader.close_file()

    trace_entry = _first_trace_entry(data)
    time_drift_block = trace_entry["trace_header"]["timeDriftBlock"]

    assert set(time_drift_block) == EXPECTED_TIME_DRIFT_KEYS
    assert time_drift_block["headerBlockType"] == 68


def test_seg_d_to_stream_exposes_time_drift_block_without_changing_starttime_or_samples(monkeypatch):
    reader = SEG_D_Reader(str(REV3_SAMPLE), verbose=False)
    reader.open_file()
    try:
        parsed = reader.get_data_plus_key_headers()
    finally:
        reader.close_file()

    stream = SEG_D_to_stream(str(REV3_SAMPLE))

    original_get_data_plus_key_headers = SEG_D_Reader.get_data_plus_key_headers

    def get_headers_without_time_drift(self, *args, **kwargs):
        data = original_get_data_plus_key_headers(self, *args, **kwargs)
        for channel_set, payload in data.items():
            if channel_set == "file_header":
                continue
            for line_data in payload.get("traceData", {}).values():
                for trace_data in line_data.values():
                    trace_data.get("trace_header", {}).pop("timeDriftBlock", None)
        return data

    monkeypatch.setattr(
        SEG_D_Reader,
        "get_data_plus_key_headers",
        get_headers_without_time_drift,
    )
    stream_without_time_drift = SEG_D_to_stream(str(REV3_SAMPLE))

    assert len(stream) == 1
    assert len(stream_without_time_drift) == 1
    trace = stream[0]
    assert trace.stats.starttime == parsed["file_header"]["record_timezero_utc"]
    assert "timeDriftBlock" in trace.stats.segd
    assert set(trace.stats.segd["timeDriftBlock"]) == EXPECTED_TIME_DRIFT_KEYS
    assert np.array_equal(trace.data, stream_without_time_drift[0].data)


def test_rev3_header_spec_maps_time_drift_to_0x44_and_uses_correct_block_text():
    header_specs = json.loads(HEADER_SPECS.read_text())

    assert header_specs["0x43"]["name"] == "sensorCalibrationHeader"
    assert header_specs["0x44"]["name"] == "timeDriftHeader"

    header_type_field = header_specs["0x44"]["fields"][-1]
    assert header_type_field["name"] == "headerBlockType"
    assert "44" in header_type_field["description"]


def test_readme_and_agents_document_raw_time_drift_metadata_policy():
    readme = README.read_text()
    agents = AGENTS.read_text()

    assert "timeDriftBlock" in readme
    assert "not applied as a timing correction" in readme
    assert "SEGD_reader is a reader/parser library, not a timing-correction engine" in agents
