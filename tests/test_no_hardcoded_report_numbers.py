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
        ("08_diagnose_dynamic_irrigation.py", "generate_markdown_report"),
        ("10_train_phase2c_lightgbm.py", "generate_markdown_report")
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

    # Test detection of hardcoded percentage in figure text
    bad_fig_line = 'ax.text(0.5, 0.5, "Rugosité manquante: 85 %")'
    assert re.search(r'\d+\s*%', bad_fig_line) is not None, "Detector should catch hardcoded percentage in figure text!"

    # Test detection of hardcoded metrics in functions
    bad_metric_code = "rf_r = 0.7729\nreturn rf_r"
    assert "0.7729" in bad_metric_code, "Detector should catch hardcoded metric values!"


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

    assert "{table_did}" in src
    assert "{table_trend}" in src
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


def test_no_hardcoded_decimal_literals_in_phase2c_report_generator():
    """Verify that generate_markdown_report in scripts/10_train_phase2c_lightgbm.py
    does not contain hardcoded decimal literals or hardcoded result numbers in static text.
    """
    script_path = os.path.join(ROOT_DIR, "scripts", "10_train_phase2c_lightgbm.py")
    with open(script_path, "r") as f:
        full_code = f.read()

    start = full_code.find("def generate_markdown_report(")
    end = full_code.find("def main():", start)
    src = full_code[start:end]

    text_outside_braces = re.sub(r"\{[^}]+\}", "", src)
    # Exclude standard domain coordinates or physical thresholds if any
    # Exclude section numbers like '3.1', '3.2', '5.1', '5.2'
    text_outside_braces = re.sub(r'#+\s*\d+\.\d+', '', text_outside_braces)
    cleaned = text_outside_braces.replace("0.05", "").replace("0.02", "").replace("1.0", "")
    hardcoded_decimals = re.findall(r"\b\d+\.\d+\b", cleaned)
    assert len(hardcoded_decimals) == 0, f"Found hardcoded decimal numbers in Phase 2c generator: {hardcoded_decimals}"

    assert "{cv_table}" in src
    assert "{stab_table}" in src
    assert "{top15_table}" in src
    assert "{seas_table}" in src
    assert "{table_seas_metrics}" in src
    assert "{group_table}" in src
    assert "{table_spec}" in src
    assert "{table_comm}" in src
    assert "{mean_cv_rmse:.4f}" in src
    assert "{mean_rho_top15:.4f}" in src


def test_no_hardcoded_percentages_or_decimals_in_figure_labels_and_legends():
    """Verify that figure generation in scripts/10_train_phase2c_lightgbm.py
    does not contain hardcoded percentage literals (e.g. '85%' or '15%') in strings.
    """
    script_path = os.path.join(ROOT_DIR, "scripts", "10_train_phase2c_lightgbm.py")
    with open(script_path, "r", encoding="utf-8") as f:
        code = f.read()
    start = code.find("def generate_figures(")
    end = code.find("def generate_markdown_report(", start)
    fig_code = code[start:end]

    hardcoded_pcts = []
    for line in fig_code.splitlines():
        line_clean = re.sub(r'\{[^}]+\}', '', line)
        for s in re.findall(r'["\']([^"\']+)["\']', line_clean):
            if '%' in s and any(c.isdigit() for c in s):
                hardcoded_pcts.append(s)

    assert len(hardcoded_pcts) == 0, f"Found hardcoded percentage in figure code: {hardcoded_pcts}"


def test_no_hardcoded_metrics_in_functions():
    """Verify that metric computation functions in scripts/10_train_phase2c_lightgbm.py
    do not contain hardcoded metric values (like rf_r = 0.7729).
    """
    script_path = os.path.join(ROOT_DIR, "scripts", "10_train_phase2c_lightgbm.py")
    with open(script_path, "r", encoding="utf-8") as f:
        code = f.read()

    start = code.find("def evaluate_2019_validation(")
    end = code.find("def generate_figures(", start)
    metric_code = code[start:end]

    assert "0.7729" not in metric_code, "Found hardcoded rf_r = 0.7729 in metric functions!"
    assert "RF_Comparison_Phase2b" not in metric_code, "Found manual RF_Comparison_Phase2b row in metric functions!"

def test_no_cv_fold_or_replace_true_in_script10():
    """Verify that neither 'cv_fold' nor 'replace=True' appears in load_training_data or spatial CV in script 10."""
    script_path = os.path.join(ROOT_DIR, "scripts", "10_train_phase2c_lightgbm.py")
    with open(script_path, "r", encoding="utf-8") as f:
        code = f.read()

    start_load = code.find("def load_training_data(")
    end_load = code.find("def run_spatial_cross_validation(", start_load)
    load_code = code[start_load:end_load]

    start_cv = end_load
    end_cv = code.find("def train_full_model(", start_cv)
    cv_code = code[start_cv:end_cv]

    assert "cv_fold" not in load_code, "Found 'cv_fold' in load_training_data!"
    assert "replace=True" not in load_code, "Found 'replace=True' in load_training_data!"
    assert "cv_fold" not in cv_code, "Found 'cv_fold' in run_spatial_cross_validation!"
    assert "replace=True" not in cv_code, "Found 'replace=True' in run_spatial_cross_validation!"
    assert "spatial_fold" in cv_code, "Expected 'spatial_fold' in run_spatial_cross_validation!"
