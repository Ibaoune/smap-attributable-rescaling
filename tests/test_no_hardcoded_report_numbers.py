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
    """Verify that generate_refined_synthesis_report in scripts/07_run_phase2a_baselines.py
    does not contain hardcoded decimal literals or hardcoded result numbers in static text.
    """
    script_path = os.path.join(ROOT_DIR, "scripts", "07_run_phase2a_baselines.py")
    with open(script_path, "r") as f:
        full_code = f.read()

    start = full_code.find("def generate_refined_synthesis_report(")
    end = full_code.find("def main():", start)
    src = full_code[start:end]

    # No hardcoded decimal numbers outside format expressions {xxx}
    text_outside_braces = re.sub(r'\{[^}]+\}', '', src)
    cleaned = text_outside_braces.replace("0.02", "").replace("1e-5", "").replace("0.05", "").replace("0.01", "").replace("1.0", "")
    hardcoded_decimals = re.findall(r'\b\d+\.\d+\b', cleaned)
    assert len(hardcoded_decimals) == 0, f"Found hardcoded decimal numbers in Phase 2a generator: {hardcoded_decimals}"

    # Verify dynamic interpolation tokens
    assert "{glob['bias_b0_spec']:+.4f}" in src
    assert "{glob['rmse_b0_spec']:.4f}" in src
    assert "{glob['rmse_b1_spec']:.4f}" in src
    assert "{glob['rmse_b2_spec']:.4f}" in src
    assert "{daobs_table}" in src
    assert "{pop_table}" in src
    assert "{contrast_table}" in src


def test_no_interpretive_phrasing_in_report_generators():
    """Verify that neither script 05 nor script 07 contains subjective/interpretive phrasing in report templates."""
    banned_phrases = [
        "quasi-parfaitement",
        "validé à 100",
        "rejet massif",
        "parfaitement",
        "fidèlement"
    ]
    for script_name, func_name in [
        ("05_assemble_modeling_tables.py", "generate_synthesis_report"),
        ("07_run_phase2a_baselines.py", "generate_refined_synthesis_report"),
        ("08_diagnose_dynamic_irrigation.py", "generate_markdown_report")
    ]:
        script_path = os.path.join(ROOT_DIR, "scripts", script_name)
        with open(script_path, "r") as f:
            code = f.read()
        start = code.find(f"def {func_name}(")
        end = code.find("def main():", start)
        func_src = code[start:end]

        for phrase in banned_phrases:
            assert phrase not in func_src.lower(), f"Found banned interpretive phrase '{phrase}' in {script_name} ({func_name})!"


def test_detector_catches_hardcoded_numbers_and_interpretations():
    """Verify that text with hardcoded decimals or interpretive phrasing is properly caught by test rules."""
    bad_source = 'rf.write("> La moyenne de d passe de -0,0362 m3/m3 à -0,0588 m3/m3, quasi-parfaitement validé à 100 %.")'
    matches = re.findall(r'rf\.write\([^)]*?\b\d+,\d+\b[^)]*?\)', bad_source)
    assert len(matches) > 0, "Detector should catch hardcoded decimal literals!"
    
    banned = ["quasi-parfaitement", "validé à 100"]
    found_banned = [p for p in banned if p in bad_source]
    assert len(found_banned) == 2, "Detector should catch banned interpretive phrases!"


def test_no_hardcoded_decimal_literals_in_phase2a_bis_report_generator():
    """Verify that generate_markdown_report in scripts/08_diagnose_dynamic_irrigation.py
    does not contain hardcoded decimal literals or hardcoded result numbers in static text.
    """
    script_path = os.path.join(ROOT_DIR, "scripts", "08_diagnose_dynamic_irrigation.py")
    with open(script_path, "r") as f:
        full_code = f.read()

    start = full_code.find("def generate_markdown_report(")
    end = full_code.find("def main():", start)
    src = full_code[start:end]

    text_outside_braces = re.sub(r'\{[^}]+\}', '', src)
    cleaned = text_outside_braces.replace("1.0", "").replace("0.03", "").replace("0.02", "").replace("0.04", "").replace("0.5", "")
    hardcoded_decimals = re.findall(r'\d+\.\d+', cleaned)
    assert len(hardcoded_decimals) == 0, f"Found hardcoded decimal numbers in Phase 2a-bis generator: {hardcoded_decimals}"

    assert "{table_cl_jja}" in src
    assert "{table_loam_jja}" in src
    assert "{table_months}" in src
    assert "{table_sens}" in src


def test_no_hardcoded_decimal_literals_in_phase2b_report_generator():
    """Verify that generate_rf_report in scripts/09_train_phase2b_random_forest.py
    does not contain hardcoded decimal literals or hardcoded result numbers in static text.
    """
    script_path = os.path.join(ROOT_DIR, "scripts", "09_train_phase2b_random_forest.py")
    with open(script_path, "r") as f:
        full_code = f.read()

    start = full_code.find("def generate_rf_report(")
    end = full_code.find("def main():", start)
    src = full_code[start:end]

    text_outside_braces = re.sub(r"\{[^}]+\}", "", src)
    cleaned = text_outside_braces.replace("1.0", "").replace("0.02", "").replace("0.01", "")
    hardcoded_decimals = re.findall(r"\b\d+\.\d+\b", cleaned)
    assert len(hardcoded_decimals) == 0, f"Found hardcoded decimal numbers in Phase 2b generator: {hardcoded_decimals}"

    assert "{cv_table}" in src
    assert "{shap_table}" in src
    assert "{pop_table}" in src
    assert "{glob_metrics_table}" in src
    assert "{season_table}" in src
    assert "{anom_table}" in src
    assert "{contrast_table}" in src
    assert "{wlt_table}" in src
