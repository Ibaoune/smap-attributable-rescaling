import re
import os
import pytest

ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))


def test_no_hardcoded_decimal_literals_in_phase1_report_generator():
    """Verify that generate_synthesis_report in scripts/05_assemble_modeling_tables.py
    does not contain hardcoded decimal literals or hardcoded result numbers in static text.
    """
    script_path = os.path.join(ROOT_DIR, "scripts", "05_assemble_modeling_tables.py")
    with open(script_path, "r") as f:
        full_code = f.read()

    start = full_code.find("def generate_synthesis_report(")
    end = full_code.find("def main():", start)
    src = full_code[start:end]

    # 1. No French decimal comma numbers (e.g., -0,0362 or 81,5)
    french_decimals = re.findall(r'rf\.write\([^)]*?\b\d+,\d+\b[^)]*?\)', src)
    assert len(french_decimals) == 0, f"Found hardcoded French decimal numbers in rf.write: {french_decimals}"

    # 2. No hardcoded MSWEP table lines
    hardcoded_mswep_rows = re.findall(r'rf\.write\([^)]*?\|\s*2016\s*\|\s*\d+\.\d+\s*\|[^)]*?\)', src)
    assert len(hardcoded_mswep_rows) == 0, f"Found hardcoded MSWEP table rows: {hardcoded_mswep_rows}"

    # 3. Dynamic GMIA note presence
    assert "{sub_0:" in src and "{sub_gt30:" in src, "GMIA note must be dynamically formatted from data variables!"

    # 4. Dynamic roughness percentage presence
    assert "{pct_7_16:" in src, "Roughness concentration percentage must be dynamically computed!"

    # 5. Dynamic QC table rule loop presence
    assert "{prec:.2f}%" in src and "{share_lis:.2f}%" in src, "QC table metrics must be dynamically computed!"


def test_no_hardcoded_decimal_literals_in_phase2a_report_generator():
    """Verify that generate_synthesis_report in scripts/07_run_phase2a_baselines.py
    does not contain hardcoded decimal literals or hardcoded result numbers in static text.
    """
    script_path = os.path.join(ROOT_DIR, "scripts", "07_run_phase2a_baselines.py")
    with open(script_path, "r") as f:
        full_code = f.read()

    start = full_code.find("def generate_synthesis_report(")
    end = full_code.find("def main():", start)
    src = full_code[start:end]

    # No hardcoded decimal numbers like -0.0265 or 0.0597 or 46.40 in template strings outside format braces
    # Remove format expressions {xxx}
    text_outside_braces = re.sub(r'\{[^}]+\}', '', src)
    # Check for hardcoded decimal literals (e.g. 0.0265, -0.0265, 46.40%)
    # Note: 0.02 in formula explanation "\theta_{wlt} + 0.02" is a definition parameter, not an empirical result.
    # Exclude the mathematical definition string "0.02"
    cleaned = text_outside_braces.replace("0.02", "").replace("1e-5", "")
    hardcoded_decimals = re.findall(r'\b\d+\.\d+\b', cleaned)
    assert len(hardcoded_decimals) == 0, f"Found hardcoded decimal numbers in Phase 2a generator: {hardcoded_decimals}"

    # Verify dynamic interpolation tokens
    assert "{glob['b0_bias']:+.4f}" in src
    assert "{glob['b0_rmse']:.4f}" in src
    assert "{glob['b1_rmse']:.4f}" in src
    assert "{glob['b2_rmse']:.4f}" in src
    assert "{daobs_table}" in src
    assert "{strata_table}" in src


def test_detector_catches_hardcoded_numbers():
    """Verify that a text with hardcoded decimals is properly caught by our detection rule."""
    bad_source = 'rf.write("> La moyenne de d passe de -0,0362 m3/m3 à -0,0588 m3/m3.")'
    matches = re.findall(r'rf\.write\([^)]*?\b\d+,\d+\b[^)]*?\)', bad_source)
    assert len(matches) > 0, "Detector should catch hardcoded decimal literals!"
