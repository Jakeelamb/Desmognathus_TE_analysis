from __future__ import annotations

import csv
import hashlib
import json
import re
import subprocess
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
NOTEBOOK_DIR = ROOT / "notebooks/research_review"
README = NOTEBOOK_DIR / "README.md"
BUNDLE_MANIFEST = NOTEBOOK_DIR / "research_review_manifest.json"
FIGURE_DIR = ROOT / "results/figures/research_review"
DATA_DIR = ROOT / "results/data/research_review"
FIGURE_MANIFEST = DATA_DIR / "audited_historical_style_figures_v1.manifest.json"
FROZEN_REGISTRY = DATA_DIR / "frozen_input_registry.csv"
SOURCE_R_RUNNER = ROOT / "scripts/processing/build_audited_historical_style_figures.R"
CORRECTED_LTR_DIR = ROOT / "results/data/corrected/ectopic_ltr30"
CORRECTED_LTR_ELEMENTS = CORRECTED_LTR_DIR / "ectopic_element_metrics_ltr30_v1.csv"
CORRECTED_LTR_EXCLUSIONS = CORRECTED_LTR_DIR / "ectopic_excluded_elements_ltr30_v1.csv"
CORRECTED_LTR_MANIFEST = CORRECTED_LTR_DIR / "ectopic_ltr30_v1.manifest.json"
LTR_ARTIFACT_SCREEN = DATA_DIR / "ltr_mapping_artifact_screen_ltr30_v1.csv"
HISTORICAL_LTR30_BASENAME = (
    "ectopic_recombination_filtered_3000bp_5+domains_no_unknown_species.csv"
)
RETRY_MANIFEST = (
    ROOT
    / "results/data/corrected/te34_replicate_averaged"
    / "te34_replicate_averaged_v1.manifest.json"
)

EXPECTED_NOTEBOOKS = (
    "01_phylogeny_tree_trimming.ipynb",
    "02_data_tables_and_provenance.ipynb",
    "03_repeat_analysis_te34.ipynb",
    "04_ltr_deletion_footprint.ipynb",
    "05_cell_modeling_and_measurement.ipynb",
    "06_genome_size_estimation.ipynb",
    "07_cell_nucleus_genome_path_analysis.ipynb",
    "08_integrated_phylogenetic_path_analysis.ipynb",
)

GENOME_DATA_DIR = DATA_DIR / "genome_size_estimation"
GENOME_FIGURE_DIR = FIGURE_DIR / "genome_size_estimation"
GENOME24_PATH_DATA_DIR = DATA_DIR / "cell_nucleus_genome_path"
GENOME24_PATH_FIGURE_DIR = FIGURE_DIR / "cell_nucleus_genome_path"

EXPECTED_FIGURE_STEMS = (
    "assembly_quality_te34_v1",
    "dnapipete_quality_te34_v1",
    "te_mass_composition_te34_v1",
    "te_mean_order_te34_v1",
    "te_mean_superfamily_te34_v1",
    "te_diversity_te34_v1",
    "te_shannon_cross_taxon_context_te34_v1",
    "te_shannon_genome_size_scatter_te18_v1",
    "te_pca_species_te34_v1",
    "te_pca_scree_te34_v1",
    "te_pca_elbow_te34_v1",
    "te_pca_silhouette_te34_v1",
    "te_pca_clusters_te34_v1",
    "te_pca_clusters_phylogeny_te34_v1",
    "ltr_terminal_internal_violin_ltr30_v1",
    "ltr_terminal_internal_artifact_screened_ltr30_v1",
    "ltr_terminal_internal_log_ltr30_v1",
)

GENOME_SHANNON_SCATTER = (
    DATA_DIR / "te_shannon_genome_size_scatter_te18_v1.csv"
)
GENOME_SHANNON_MATCH_AUDIT = (
    DATA_DIR / "te_shannon_genome_size_match_audit_te34_v1.csv"
)
GENOME_SHANNON_STATS = DATA_DIR / "te_shannon_genome_size_stats_te18_v1.txt"
SUPERFAMILY_PREVALENCE = DATA_DIR / "te_superfamily_prevalence_te34_v1.csv"
SUPERFAMILY_PCA_MATRIX = DATA_DIR / "te_pca_superfamily_feature_matrix_te34_v1.csv"

PANEL_FILES = {
    "TE34": (
        ROOT / "path_analysis/data/derived/panels/study_te_resource_panel34_v1.csv",
        34,
    ),
    "Cell21": (
        ROOT / "path_analysis/data/derived/panels/study_cell_linked_panel21_v1.csv",
        21,
    ),
    "path18": (
        ROOT / "path_analysis/data/derived/panels/study_integrated_path_panel18_v1.csv",
        18,
    ),
}

EXPENSIVE_TOOL_NAMES = (
    "repeatmasker",
    "repeatmodeler",
    "dnapipete",
    "cellpose",
    "tesorter",
    "fastqc",
    "multiqc",
    "snakemake",
    "nextflow",
)


def require_file(test: unittest.TestCase, path: Path) -> None:
    """Give a useful assertion in the intentional pre-build red state."""

    test.assertTrue(path.is_file(), f"required file is missing: {path.relative_to(ROOT)}")


def read_notebook(test: unittest.TestCase, name: str) -> dict[str, object]:
    path = NOTEBOOK_DIR / name
    require_file(test, path)
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as error:
        test.fail(f"invalid notebook JSON in {name}: {error}")
    raise AssertionError("unreachable")


def cell_source(cell: dict[str, object]) -> str:
    source = cell.get("source", "")
    if isinstance(source, list):
        return "".join(str(part) for part in source)
    return str(source)


def notebook_source(notebook: dict[str, object], cell_type: str | None = None) -> str:
    cells = notebook.get("cells", [])
    if not isinstance(cells, list):
        return ""
    return "\n".join(
        cell_source(cell)
        for cell in cells
        if isinstance(cell, dict)
        and (cell_type is None or cell.get("cell_type") == cell_type)
    )


def notebook_output_text(notebook: dict[str, object]) -> str:
    chunks: list[str] = []
    cells = notebook.get("cells", [])
    if not isinstance(cells, list):
        return ""
    for cell in cells:
        if not isinstance(cell, dict):
            continue
        outputs = cell.get("outputs", [])
        if not isinstance(outputs, list):
            continue
        for output in outputs:
            if not isinstance(output, dict):
                continue
            text = output.get("text")
            if isinstance(text, list):
                chunks.extend(str(part) for part in text)
            elif text is not None:
                chunks.append(str(text))
            data = output.get("data", {})
            if isinstance(data, dict):
                for mime in ("text/plain", "text/html", "text/markdown"):
                    value = data.get(mime)
                    if isinstance(value, list):
                        chunks.extend(str(part) for part in value)
                    elif value is not None:
                        chunks.append(str(value))
    return "\n".join(chunks)


def normalized_bool(value: object) -> bool:
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)):
        return value != 0
    return str(value).strip().lower() in {"true", "t", "yes", "y", "1"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def species_values(path: Path) -> list[str]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames or []
        species_field = next(
            (
                field
                for field in fields
                if re.sub(r"[^a-z]", "", field.lower())
                in {"species", "speciesname", "taxon", "taxonname"}
            ),
            None,
        )
        if species_field is None:
            return []
        return [str(row.get(species_field, "")).strip().lower() for row in reader]


def expensive_invocations(code: str) -> list[str]:
    """Detect execution surfaces while allowing tool names in frozen data paths."""

    findings: list[str] = []
    process_patterns = (
        r"\bsubprocess\s*\.\s*(?:run|popen|call|check_call|check_output)\s*\(",
        r"\bos\s*\.\s*(?:system|popen)\s*\(",
        r"\bget_ipython\s*\(\s*\)\s*\.\s*system\s*\(",
        r"\b(?:system|popen)\s*\(\s*['\"](?:rscript|repeatmasker|repeatmodeler|dnapipete)",
    )
    for pattern in process_patterns:
        if re.search(pattern, code, flags=re.IGNORECASE):
            findings.append(pattern)

    tool_alternation = "|".join(re.escape(name) for name in EXPENSIVE_TOOL_NAMES)
    shell_pattern = rf"(?im)^\s*(?:!|%%?bash\b|%%?sh\b).*\b(?:{tool_alternation})\b"
    for match in re.finditer(shell_pattern, code):
        findings.append(match.group(0).strip())

    direct_api_patterns = (
        r"(?im)^\s*(?:from\s+cellpose\b|import\s+cellpose\b)",
        r"(?im)^\s*(?:from\s+snakemake\b|import\s+snakemake\b)",
    )
    for pattern in direct_api_patterns:
        if re.search(pattern, code):
            findings.append(pattern)
    return findings


class ResearchReviewBundleTests(unittest.TestCase):
    maxDiff = None

    def test_exactly_eight_numbered_domain_notebooks_exist(self) -> None:
        self.assertTrue(
            NOTEBOOK_DIR.is_dir(),
            f"required directory is missing: {NOTEBOOK_DIR.relative_to(ROOT)}",
        )
        actual = sorted(path.name for path in NOTEBOOK_DIR.glob("*.ipynb"))
        self.assertEqual(actual, sorted(EXPECTED_NOTEBOOKS))
        self.assertTrue(all(re.match(r"^0[1-8]_", name) for name in actual))

    def test_superseded_presentation_tree_is_removed(self) -> None:
        presentation = ROOT / "notebooks/presentation"
        self.assertFalse(
            presentation.exists(),
            "superseded presentation notebooks must not coexist with the canonical bundle",
        )

    def test_genome_size_estimation_is_a_complete_standalone_notebook(self) -> None:
        notebook = read_notebook(self, "06_genome_size_estimation.ipynb")
        source = notebook_source(notebook)
        self.assertIn("Finalized Quality-Matched Genome-Size Analysis", source)
        self.assertIn("build_frozen_genome_iod_notebook.py", source)
        self.assertIn("build_frozen_genome_iod_phylogeny_figure.py", source)
        self.assertIn("D. fuscus", source)
        self.assertIn("conditional", source.lower())
        self.assertRegex(source, r"not an independent\s+direct C-value")

        for name in (
            "species_relative_genome_iod_summary.csv",
            "image_relative_genome_iod_summary.csv",
            "frozen_quality_balance.csv",
            "iod_quality_residual_diagnostics.csv",
            "phylogeny_genome_nucleus_cell_summary.csv",
            "phylogeny_genome_nucleus_cell_correlations.csv",
            "analysis_manifest.json",
            "phylogeny_genome_nucleus_cell_manifest.json",
        ):
            require_file(self, GENOME_DATA_DIR / name)

        for name in (
            "01_all_raw_nuclear_iod.png",
            "02_relative_iod_estimates.png",
            "03_image_level_iod.png",
            "04_iod_components.png",
            "05_aggregation_sensitivity.png",
            "06_quality_diagnostics.png",
            "07_measured_phylogeny_genome_nucleus_cell.png",
            "07_measured_phylogeny_genome_nucleus_cell.pdf",
            "08_pairwise_genome_nucleus_cell_relationships.png",
            "08_pairwise_genome_nucleus_cell_relationships.pdf",
        ):
            require_file(self, GENOME_FIGURE_DIR / name)

        self.assertFalse((NOTEBOOK_DIR / "06_genome_size_estimation.executed.ipynb").exists())
        self.assertFalse((NOTEBOOK_DIR / "06_genome_size_estimation.html").exists())

    def test_genome24_path_analysis_is_full_and_does_not_force_a_causal_winner(self) -> None:
        notebook = read_notebook(self, "07_cell_nucleus_genome_path_analysis.ipynb")
        source = notebook_source(notebook)
        self.assertIn("Finalized 24-species", source)
        self.assertIn("25", source)
        self.assertIn("11", source)
        self.assertIn("Markov-equivalent", source)
        self.assertIn("250", source)
        self.assertIn("200", source)
        self.assertIn("100", source)
        self.assertRegex(source, r"(?i)do not identify a unique causal pathway")

        analysis = json.loads(
            (GENOME24_PATH_DATA_DIR / "analysis_manifest.json").read_text()
        )
        self.assertEqual(analysis["run_scope"], "full_release")
        self.assertEqual(analysis["measurement_bootstrap_replicates"], 250)
        self.assertEqual(analysis["published_bootstrap_trees_analyzed"], 200)
        self.assertEqual(analysis["simulation_replicates_per_class_and_regime"], 100)
        self.assertEqual(analysis["failure_count"], 0)
        self.assertEqual(analysis["primary_species"], 24)
        self.assertEqual(analysis["primary_top_model"], "cell_bridge")
        self.assertEqual(analysis["primary_second_model"], "cell_collider")
        self.assertGreater(analysis["primary_top_two_delta_cicc"], 2.0)
        self.assertFalse(analysis["causal_direction_identified"])

        figure_manifest = json.loads(
            (
                GENOME24_PATH_DATA_DIR
                / "cell_nucleus_genome_path_figure_manifest.json"
            ).read_text()
        )
        self.assertTrue(figure_manifest["release_gates_passed"])
        self.assertEqual(len(figure_manifest["figures"]), 8)
        for record in figure_manifest["figures"]:
            for key in ("png", "pdf"):
                path = ROOT / record[key]
                require_file(self, path)
                self.assertEqual(record[f"{key}_sha256"], sha256(path))

    def test_bundle_manifest_hashes_the_exact_canonical_notebooks(self) -> None:
        require_file(self, BUNDLE_MANIFEST)
        manifest = json.loads(BUNDLE_MANIFEST.read_text(encoding="utf-8"))
        self.assertEqual(manifest["status"], "pass")
        self.assertEqual(set(manifest["notebooks"]), set(EXPECTED_NOTEBOOKS))
        self.assertTrue(manifest["gates"]["all_code_cells_executed"])
        self.assertTrue(manifest["gates"]["no_error_outputs"])
        self.assertTrue(
            manifest["gates"]["no_expensive_or_external_execution_detected_in_code"]
        )
        for name, record in manifest["notebooks"].items():
            self.assertEqual(record["sha256"], sha256(NOTEBOOK_DIR / name), name)

        manifest_figures = {
            record["path"]: record["sha256"] for record in manifest["figure_outputs"]
        }
        actual_figures = {
            str(path.relative_to(ROOT))
            for path in FIGURE_DIR.rglob("*")
            if path.is_file() and path.suffix.lower() in {".png", ".pdf", ".svg"}
        }
        self.assertEqual(set(manifest_figures), actual_figures)
        for relative, expected_hash in manifest_figures.items():
            self.assertEqual(expected_hash, sha256(ROOT / relative), relative)

    def test_notebooks_are_fully_executed_and_have_no_error_outputs(self) -> None:
        for name in EXPECTED_NOTEBOOKS:
            with self.subTest(notebook=name):
                notebook = read_notebook(self, name)
                cells = notebook.get("cells", [])
                self.assertIsInstance(cells, list, name)
                code_cells = [
                    cell
                    for cell in cells
                    if isinstance(cell, dict) and cell.get("cell_type") == "code"
                ]
                self.assertGreaterEqual(len(code_cells), 1, f"{name} has no code cells")
                for index, cell in enumerate(code_cells, start=1):
                    execution_count = cell.get("execution_count")
                    self.assertIsInstance(
                        execution_count,
                        int,
                        f"{name} code cell {index} was not executed",
                    )
                    self.assertGreater(
                        execution_count,
                        0,
                        f"{name} code cell {index} has an invalid execution count",
                    )
                    outputs = cell.get("outputs", [])
                    self.assertIsInstance(outputs, list, name)
                    errors = [
                        output
                        for output in outputs
                        if isinstance(output, dict) and output.get("output_type") == "error"
                    ]
                    self.assertEqual(errors, [], f"{name} code cell {index}")
                self.assertNotIn(
                    "Traceback (most recent call last)",
                    notebook_output_text(notebook),
                    name,
                )

    def test_notebook_cells_do_not_launch_expensive_upstream_tools(self) -> None:
        for name in EXPECTED_NOTEBOOKS:
            with self.subTest(notebook=name):
                notebook = read_notebook(self, name)
                code = notebook_source(notebook, "code")
                self.assertEqual(expensive_invocations(code), [], name)

    def test_frozen_registry_hashes_every_portable_dependency(self) -> None:
        require_file(self, FROZEN_REGISTRY)
        rows = csv_rows(FROZEN_REGISTRY)
        self.assertGreater(len(rows), 0)
        required_columns = {
            "artifact_kind",
            "source_path",
            "frozen_path",
            "size_bytes",
            "source_sha256",
            "frozen_sha256",
            "copied_at_utc",
            "source_present_at_build",
            "source_git_status",
        }
        self.assertTrue(required_columns.issubset(rows[0]), "frozen registry schema")
        self.assertEqual(len({row["source_path"] for row in rows}), len(rows))
        self.assertEqual(len({row["frozen_path"] for row in rows}), len(rows))
        self.assertEqual(
            {row["artifact_kind"] for row in rows},
            {"data_or_tree", "figure"},
        )

        for row in rows:
            source_relative = Path(row["source_path"])
            frozen_relative = Path(row["frozen_path"])
            with self.subTest(snapshot=row["frozen_path"]):
                self.assertFalse(source_relative.is_absolute())
                self.assertFalse(frozen_relative.is_absolute())
                self.assertNotIn("..", source_relative.parts)
                self.assertNotIn("..", frozen_relative.parts)
                frozen = ROOT / frozen_relative
                require_file(self, frozen)
                self.assertEqual(int(row["size_bytes"]), frozen.stat().st_size)
                self.assertRegex(row["source_sha256"], r"^[0-9a-f]{64}$")
                self.assertEqual(row["source_sha256"], row["frozen_sha256"])
                self.assertEqual(row["frozen_sha256"], sha256(frozen))
                self.assertTrue(row["copied_at_utc"])
                self.assertIn(
                    row["source_git_status"],
                    {"tracked", "ignored", "untracked", "source_absent"},
                )
                if normalized_bool(row["source_present_at_build"]):
                    source = ROOT / source_relative
                    require_file(self, source)
                    self.assertEqual(row["source_sha256"], sha256(source))
                if row["artifact_kind"] == "data_or_tree":
                    self.assertTrue(
                        str(frozen_relative).startswith(
                            "results/data/research_review/frozen_inputs/"
                        )
                    )
                else:
                    self.assertTrue(
                        str(frozen_relative).startswith("results/figures/research_review/")
                    )
                ignored = subprocess.run(
                    ["git", "check-ignore", "-q", "--", str(frozen_relative)],
                    cwd=ROOT,
                    check=False,
                ).returncode == 0
                self.assertFalse(ignored, f"portable snapshot is ignored: {frozen_relative}")

    def test_notebook_file_dependencies_resolve_through_frozen_snapshots(self) -> None:
        registry = csv_rows(FROZEN_REGISTRY)
        source_to_frozen = {row["source_path"]: row["frozen_path"] for row in registry}
        file_suffixes = {".csv", ".json", ".nwk", ".nex", ".tre", ".txt", ".png", ".pdf"}
        for name in EXPECTED_NOTEBOOKS:
            with self.subTest(notebook=name):
                notebook = read_notebook(self, name)
                code = notebook_source(notebook, "code")
                for required_token in (
                    "frozen_input_registry.csv",
                    "FROZEN_ARTIFACTS",
                    "resolve_artifact",
                ):
                    self.assertIn(required_token, code, f"{name}: {required_token}")
                literal_paths = set(
                    re.findall(
                        r'''["']((?:input_data|path_analysis/data|results)/[^"']+)["']''',
                        code,
                    )
                )
                upstream_files = {
                    path
                    for path in literal_paths
                    if Path(path).suffix.lower() in file_suffixes
                    and not path.startswith("results/data/research_review/")
                    and not path.startswith("results/figures/research_review/")
                }
                missing = sorted(upstream_files - set(source_to_frozen))
                self.assertEqual(missing, [], f"{name} has unsnapshotted dependencies")
                for source in upstream_files:
                    self.assertTrue((ROOT / source_to_frozen[source]).is_file(), source)

    def test_saved_notebooks_have_no_live_localhost_viewer_url(self) -> None:
        live_url = re.compile(r"https?://(?:127\.0\.0\.1|localhost):\d+", re.IGNORECASE)
        for name in EXPECTED_NOTEBOOKS:
            with self.subTest(notebook=name):
                notebook = read_notebook(self, name)
                self.assertIsNone(live_url.search(notebook_output_text(notebook)), name)

        cell_notebook = read_notebook(self, "05_cell_modeling_and_measurement.ipynb")
        cell_code = notebook_source(cell_notebook, "code")
        self.assertIn("def start_optional_local_viewers", cell_code)
        self.assertIsNone(
            re.search(r"(?m)^start_optional_local_viewers\s*\(", cell_code),
            "batch execution must not start or save a live localhost viewer",
        )

    def test_readme_and_panel_tables_use_te34_cell21_and_path18(self) -> None:
        require_file(self, README)
        readme = README.read_text(encoding="utf-8")
        denominator_patterns = {
            "TE34": r"(?is)\bTE34\b.{0,100}(?:\bn\s*=\s*34\b|\b34\s+species\b)",
            "Cell21": r"(?is)\bCell21\b.{0,100}(?:\bn\s*=\s*21\b|\b21\s+species\b)",
            "path18": r"(?is)\bpath18\b.{0,100}(?:\bn\s*=\s*18\b|\b18\s+species\b)",
        }
        for panel, pattern in denominator_patterns.items():
            self.assertRegex(readme, pattern, f"README denominator for {panel}")

        for panel, (path, expected_n) in PANEL_FILES.items():
            with self.subTest(panel=panel):
                require_file(self, path)
                species = species_values(path)
                self.assertEqual(len(species), expected_n, panel)
                self.assertEqual(len(set(species)), expected_n, f"duplicate {panel} species")
                self.assertNotIn("planiceps", species, panel)

    def test_planiceps_is_not_present_as_a_current_data_row(self) -> None:
        if DATA_DIR.is_dir():
            for path in sorted(DATA_DIR.glob("*.csv")):
                with self.subTest(table=path.name):
                    self.assertNotIn("planiceps", species_values(path))

        row_patterns = (
            # The NCBI-reported public label "Desmognathus planiceps" is
            # intentionally visible in provenance columns for the accession
            # expert-reidentified as fuscus.  Reject only a canonical/display
            # species row named planiceps.
            r"(?i)<td>\s*(?:d\.\s*)?planiceps\s*</td>",
            r"(?im)^\s*(?:d\.\s*)?planiceps\s{2,}",
        )
        for name in EXPECTED_NOTEBOOKS:
            notebook = read_notebook(self, name)
            output_text = notebook_output_text(notebook)
            for pattern in row_patterns:
                self.assertIsNone(re.search(pattern, output_text), name)

    def test_retry_runs_are_equal_weighted_within_orestes_not_double_counted(self) -> None:
        require_file(self, RETRY_MANIFEST)
        manifest = json.loads(RETRY_MANIFEST.read_text(encoding="utf-8"))
        self.assertEqual(
            manifest.get("dnapipete_run_aggregation"),
            "equal_weight_mean_of_run_level_estimates",
        )
        self.assertEqual(manifest.get("species_with_multiple_dnapipete_runs"), ["orestes"])
        self.assertIs(manifest.get("retry_runs_counted_as_additional_species"), False)

        repeat_notebook = read_notebook(self, "03_repeat_analysis_te34.ipynb")
        source = notebook_source(repeat_notebook).lower()
        self.assertIn("equal-weight", source)
        self.assertIn("within-species", source)
        self.assertIn("srx19952890", source)
        self.assertIn("srx19952890r2", source)

        require_file(self, FIGURE_MANIFEST)
        figure_manifest_text = FIGURE_MANIFEST.read_text(encoding="utf-8").lower()
        self.assertIn("equal_weight_mean_of_run_level_estimates", figure_manifest_text)
        self.assertIn("srx19952890", figure_manifest_text)
        self.assertIn("srx19952890r2", figure_manifest_text)

    def test_primary_figures_are_current_regenerations_not_historical_pngs(self) -> None:
        require_file(self, README)
        readme = README.read_text(encoding="utf-8")
        self.assertRegex(
            readme,
            r"(?is)(?:regenerat|recreat).{0,120}audited\s+(?:and\s+clean\s+)?current\s+data",
        )
        self.assertRegex(
            readme,
            r"(?is)no\s+historical\s+(?:pdf\s*/\s*png|png|snapshot).{0,80}current\s+result",
        )

        combined_source = readme + "\n"
        for name in EXPECTED_NOTEBOOKS:
            combined_source += notebook_source(read_notebook(self, name)) + "\n"

        historical_image_patterns = (
            r"(?i)(?:historical_reference|/tmp/desmo_history)[^\s'\"\)]*\.png",
            r"(?i)(?:talk_iv_restored|legacy_descriptive)[^\s'\"\)]*\.png",
            r"(?i)pdfimages\s+-(?:png|j)",
        )
        for pattern in historical_image_patterns:
            self.assertIsNone(re.search(pattern, combined_source), pattern)

        for stem in EXPECTED_FIGURE_STEMS:
            expected_reference = f"results/figures/research_review/{stem}.png"
            self.assertIn(expected_reference, combined_source, expected_reference)

    def test_relative_nuclear_iod_is_not_presented_as_absolute_pg(self) -> None:
        cell_notebook = read_notebook(self, "05_cell_modeling_and_measurement.ipynb")
        prose = (
            README.read_text(encoding="utf-8")
            + "\n"
            + notebook_source(cell_notebook, "markdown")
        )
        self.assertRegex(prose, r"(?i)relative\s+nuclear[- ]IOD")
        self.assertRegex(
            prose,
            r"(?is)(?:relative\s+nuclear[- ]IOD).{0,180}(?:not\s+(?:an?\s+)?absolute|not\s+calibrated|cannot\s+be\s+converted).{0,80}\bpg\b|"
            r"\bpg\b.{0,80}(?:not\s+(?:an?\s+)?absolute|not\s+calibrated|cannot\s+be\s+converted).{0,180}(?:relative\s+nuclear[- ]IOD)",
        )

        suspicious: list[str] = []
        for sentence in re.split(r"[.!?\n]+", prose):
            lower = sentence.lower()
            if "iod" not in lower or not re.search(r"\bpg\b", lower):
                continue
            if not re.search(
                r"\b(?:not|cannot|can't|uncalibrated|without|relative|proxy)\b",
                lower,
            ):
                suspicious.append(sentence.strip())
        self.assertEqual(
            suspicious,
            [],
            "IOD/pg statements must explicitly preserve the relative, uncalibrated boundary",
        )

    def test_genome_shannon_scatter_uses_variable_estimates_and_exact_overlap(self) -> None:
        require_file(self, GENOME_SHANNON_SCATTER)
        require_file(self, GENOME_SHANNON_MATCH_AUDIT)
        require_file(self, GENOME_SHANNON_STATS)

        scatter = csv_rows(GENOME_SHANNON_SCATTER)
        desmo = [row for row in scatter if row["group"] == "Desmognathus"]
        salamanders = [row for row in scatter if row["group"] == "Salamanders"]
        self.assertEqual(len(scatter), 28)
        self.assertEqual(len(desmo), 18)
        self.assertEqual(len(salamanders), 10)
        self.assertEqual(len({row["canonical_species"] for row in desmo}), 18)
        self.assertNotIn("planiceps", {row["canonical_species"] for row in desmo})

        genome_sizes = [float(row["genome_size_gb"]) for row in desmo]
        self.assertGreater(len(set(genome_sizes)), 10)
        self.assertGreater(max(genome_sizes) - min(genome_sizes), 5.0)
        self.assertFalse(all(abs(value - 16.0) < 1e-12 for value in genome_sizes))
        self.assertTrue(all(float(row["genome_size_ci_low_gb"]) < float(row["genome_size_gb"]) for row in desmo))
        self.assertTrue(all(float(row["genome_size_ci_high_gb"]) > float(row["genome_size_gb"]) for row in desmo))
        self.assertTrue(
            all(
                row["genome_measurement_kind"]
                == "fuscus_anchored_image_iod_genome_size_estimate"
                for row in desmo
            )
        )
        self.assertTrue(all(not normalized_bool(row["is_direct_c_value"]) for row in desmo))
        self.assertTrue(all(row["genome_support_tier"] in {"low", "limited", "medium"} for row in desmo))

        audit = csv_rows(GENOME_SHANNON_MATCH_AUDIT)
        included = [row for row in audit if normalized_bool(row["included_in_scatter"])]
        self.assertEqual(len(audit), 37)
        self.assertEqual(len(included), 18)
        self.assertEqual(
            {row["canonical_species"] for row in included},
            {row["canonical_species"] for row in desmo},
        )
        path18_species = set(species_values(PANEL_FILES["path18"][0]))
        self.assertEqual({row["canonical_species"] for row in included}, path18_species)
        self.assertTrue(all(normalized_bool(row["has_te_shannon"]) for row in included))
        self.assertTrue(all(normalized_bool(row["has_genome_size_estimate"]) for row in included))

        stats = GENOME_SHANNON_STATS.read_text(encoding="utf-8").lower()
        self.assertIn("not direct c-values", stats)
        self.assertIn("brownian-motion pgls", stats)
        self.assertRegex(stats, r"historical salamander.{0,120}excluded from all tests")
        self.assertIn("do not propagate genome-estimate or shannon-estimate uncertainty", stats)
        self.assertIn("tree uncertainty is not propagated", stats)

        manifest = json.loads(FIGURE_MANIFEST.read_text(encoding="utf-8"))
        contract = manifest["cross_taxon_shannon"]["genome_size_scatter"]
        self.assertEqual(contract["n_current_desmognathus"], 18)
        self.assertEqual(contract["n_historical_salamanders"], 10)
        self.assertIs(contract["absolute_c_value_claimed"], False)
        self.assertEqual(contract["calibration_reference_species"], "D. fuscus")
        self.assertAlmostEqual(float(contract["calibration_reference_pg"]), 16.36)
        self.assertIn("image-IOD", contract["measurement_boundary"])
        self.assertIn("not direct C-values", contract["measurement_boundary"])

        repeat_notebook = read_notebook(self, "03_repeat_analysis_te34.ipynb")
        repeat_source = notebook_source(repeat_notebook)
        self.assertIn("te_shannon_genome_size_scatter_te18_v1.csv", repeat_source)
        self.assertIn("te_shannon_genome_size_scatter_te18_v1.png", repeat_source)
        self.assertRegex(repeat_source, r"(?i)fuscus-anchored.{0,100}image-IOD")
        self.assertRegex(repeat_source, r"(?i)not direct C-values|not a direct C-value")

        frozen_sources = {row["source_path"] for row in csv_rows(FROZEN_REGISTRY)}
        self.assertIn(
            "path_analysis/data/external/derived/cellprofiler_final_species_results.csv",
            frozen_sources,
        )
        self.assertIn(
            "path_analysis/data/external/derived/cellprofiler_genome_state_summary.csv",
            frozen_sources,
        )

    def test_cross_taxon_shannon_boxplot_uses_latest_audited_te34_values(self) -> None:
        cross_path = DATA_DIR / "te_shannon_cross_taxon_context_te34_v1.csv"
        top10_path = DATA_DIR / "te_top10_superfamily_contributions_te34_v1.csv"
        require_file(self, cross_path)
        require_file(self, top10_path)
        cross = csv_rows(cross_path)
        top10 = csv_rows(top10_path)

        expected_counts = {
            "Birds": 14,
            "Bony Fishes": 13,
            "Caecilian": 4,
            "Cartilaginous fish": 3,
            "Desmognathus": 34,
            "Frogs": 10,
            "Lobe-finned fish": 3,
            "Mammals": 14,
            "Non-avian reptiles": 16,
            "Salamanders": 10,
        }
        observed_counts = {
            group: sum(row["group"] == group for row in cross)
            for group in {row["group"] for row in cross}
        }
        self.assertEqual(len(cross), 121)
        self.assertEqual(observed_counts, expected_counts)
        self.assertEqual(len(top10), 340)
        self.assertEqual(len({row["species"] for row in top10}), 34)

        contribution_sums: dict[str, float] = {}
        shannon_sums: dict[str, float] = {}
        for row in top10:
            species = row["species"]
            contribution_sums[species] = contribution_sums.get(species, 0.0) + float(
                row["top10_renormalized_proportion"]
            )
            shannon_sums[species] = shannon_sums.get(species, 0.0) + float(
                row["shannon_component"]
            )
        self.assertTrue(all(abs(value - 1.0) < 1e-10 for value in contribution_sums.values()))

        desmo = [row for row in cross if row["group"] == "Desmognathus"]
        observed_shannon = {
            row["species"][len("Desmognathus ") :]
            if row["species"].startswith("Desmognathus ")
            else row["species"]: float(
                row["shannon_top10_superfamilies"]
            )
            for row in desmo
        }
        self.assertEqual(set(observed_shannon), set(shannon_sums))
        for species, expected in shannon_sums.items():
            self.assertAlmostEqual(observed_shannon[species], expected, places=12)
        self.assertTrue(
            all("Current audited retry-averaged TE34" in row["source"] for row in desmo)
        )

        manifest = json.loads(FIGURE_MANIFEST.read_text(encoding="utf-8"))
        scope = manifest["scopes"]["cross_taxon_shannon"]
        self.assertEqual(scope["n_current_desmognathus"], 34)
        self.assertEqual(scope["n_total"], 121)
        self.assertEqual(
            manifest["dnapipete_retry_policy"]["method_id"],
            "equal_weight_mean_of_run_level_estimates",
        )

    def test_superfamily_pca_keeps_chapaev_and_dada_but_excludes_ginger(self) -> None:
        require_file(self, SUPERFAMILY_PREVALENCE)
        require_file(self, SUPERFAMILY_PCA_MATRIX)
        prevalence = {
            row["superfamily"].strip().lower(): row
            for row in csv_rows(SUPERFAMILY_PREVALENCE)
        }
        for superfamily in ("chapaev", "dada", "ginger"):
            self.assertIn(superfamily, prevalence)
            self.assertEqual(int(prevalence[superfamily]["panel_species"]), 34)

        self.assertEqual(prevalence["chapaev"]["pca_decision"], "retain")
        self.assertEqual(prevalence["dada"]["pca_decision"], "retain")
        self.assertEqual(prevalence["ginger"]["pca_decision"], "exclude")
        self.assertEqual(
            prevalence["ginger"]["decision_reason"],
            "user_declared_sparse_superfamily",
        )

        matrix = csv_rows(SUPERFAMILY_PCA_MATRIX)
        self.assertEqual(len(matrix), 34)
        fields = set(matrix[0])
        self.assertIn("Chapaev", fields)
        self.assertIn("Dada", fields)
        self.assertNotIn("Ginger", fields)

        manifest = json.loads(FIGURE_MANIFEST.read_text(encoding="utf-8"))
        contract = manifest["superfamily_pca_filter"]
        self.assertEqual(contract["panel"], "TE34")
        self.assertEqual(contract["retained_by_explicit_policy"], ["Chapaev", "Dada"])
        self.assertEqual(contract["excluded_by_explicit_policy"], ["Ginger"])

    def test_historical_style_assets_manifest_and_r_provenance_exist(self) -> None:
        require_file(self, SOURCE_R_RUNNER)
        runner_source = SOURCE_R_RUNNER.read_text(encoding="utf-8").lower()
        for provenance_token in (
            "1ad5e8c",
            "te_pca_analysis.r",
            "scripts/visualization/te_landscape_plots.r",
            "scripts/visualization/ectopic_recomb_plots.r",
        ):
            self.assertIn(provenance_token, runner_source, provenance_token)

        require_file(self, FIGURE_MANIFEST)
        manifest = json.loads(FIGURE_MANIFEST.read_text(encoding="utf-8"))
        manifest_text = json.dumps(manifest, sort_keys=True).lower()
        self.assertIs(manifest.get("expensive_upstream_tools_executed"), False)
        self.assertIn(
            "scripts/processing/build_audited_historical_style_figures.r",
            manifest_text,
        )
        self.assertRegex(manifest_text, r"audited.{0,80}current|current.{0,80}audited")

        all_notebook_source = README.read_text(encoding="utf-8") + "\n"
        for name in EXPECTED_NOTEBOOKS:
            all_notebook_source += notebook_source(read_notebook(self, name)) + "\n"
        self.assertIn(
            "scripts/processing/build_audited_historical_style_figures.R",
            all_notebook_source,
        )

        for stem in EXPECTED_FIGURE_STEMS:
            for suffix, minimum_size in ((".png", 5_000), (".pdf", 1_000)):
                path = FIGURE_DIR / f"{stem}{suffix}"
                with self.subTest(asset=path.name):
                    require_file(self, path)
                    self.assertGreater(path.stat().st_size, minimum_size, path.name)
                    expected_path = f"results/figures/research_review/{path.name}".lower()
                    self.assertIn(expected_path, manifest_text, expected_path)

    def test_manifest_and_registry_completely_cover_review_figures(self) -> None:
        manifest = json.loads(FIGURE_MANIFEST.read_text(encoding="utf-8"))
        declared_hashes: dict[str, str] = {}
        for record in manifest.get("figures", []):
            if "path" in record:
                declared_hashes[record["path"]] = record["sha256"]
                continue
            for suffix in ("png", "pdf"):
                self.assertIn(suffix, record)
                self.assertIn(f"{suffix}_sha256", record)
                declared_hashes[record[suffix]] = record[f"{suffix}_sha256"]

        registry = csv_rows(FROZEN_REGISTRY)
        for row in registry:
            if row["artifact_kind"] != "figure":
                continue
            path = row["frozen_path"]
            self.assertNotIn(path, declared_hashes, f"duplicate figure authority: {path}")
            declared_hashes[path] = row["frozen_sha256"]

        actual = {
            str(path.relative_to(ROOT))
            for path in FIGURE_DIR.iterdir()
            if path.is_file() and path.suffix.lower() in {".png", ".pdf"}
        }
        self.assertEqual(set(declared_hashes), actual)
        for relative, expected_hash in declared_hashes.items():
            with self.subTest(figure=relative):
                path = ROOT / relative
                require_file(self, path)
                self.assertRegex(expected_hash, r"^[0-9a-f]{64}$")
                self.assertEqual(sha256(path), expected_hash)

    def test_ltr30_uses_corrected_zero_aware_elements_and_audited_exclusions(self) -> None:
        require_file(self, CORRECTED_LTR_ELEMENTS)
        require_file(self, CORRECTED_LTR_EXCLUSIONS)
        require_file(self, CORRECTED_LTR_MANIFEST)
        elements = csv_rows(CORRECTED_LTR_ELEMENTS)
        exclusions = csv_rows(CORRECTED_LTR_EXCLUSIONS)
        manifest = json.loads(CORRECTED_LTR_MANIFEST.read_text(encoding="utf-8"))

        self.assertEqual(len(elements), 1086)
        self.assertEqual(len({row["species"] for row in elements}), 30)
        self.assertNotIn("planiceps", {row["species"].lower() for row in elements})
        required_element_columns = {
            "element_id",
            "species",
            "tesorter_annotation_found",
            "depth_positions_missing",
            "depth_positions_explicit_zero",
            "source_depth_status",
            "ratio_terminal_internal_all_positions",
            "ratio_terminal_internal_nonzero_only",
        }
        self.assertTrue(required_element_columns.issubset(elements[0]))
        self.assertTrue(all(row["source_depth_status"] == "usable" for row in elements))
        self.assertTrue(all(int(row["depth_positions_missing"]) == 0 for row in elements))
        self.assertTrue(all(normalized_bool(row["tesorter_annotation_found"]) for row in elements))
        self.assertGreater(
            sum(int(row["depth_positions_explicit_zero"]) for row in elements),
            0,
            "the zero-aware audit must retain explicit zero-depth positions",
        )
        self.assertTrue(
            all(float(row["ratio_terminal_internal_all_positions"]) > 0 for row in elements)
        )

        expected_exclusions = {
            "JAUEJG010675316.1_De_4598_10648": (6051, 5945),
            "JASANM010280042.1_De_179_6669": (6491, 6396),
        }
        self.assertEqual(len(exclusions), 2)
        self.assertEqual({row["element_id"] for row in exclusions}, set(expected_exclusions))
        for row in exclusions:
            expected_position, reported_position = expected_exclusions[row["element_id"]]
            self.assertEqual(row["exclusion_reason"], "truncated_depth_file_partial_final_line")
            self.assertEqual(int(row["depth_file_size_bytes"]), 245760)
            self.assertEqual(int(row["expected_terminal_position"]), expected_position)
            self.assertEqual(int(row["reported_terminal_position"]), reported_position)
            self.assertTrue(normalized_bool(row["partial_final_line"]))
            self.assertFalse(normalized_bool(row["input_modified"]))

        expected_manifest_values = {
            "analysis_scope": "audited_te_resource_panel34_ltr30",
            "n_panel_species": 34,
            "n_tabout_resources": 31,
            "n_species_with_selected_elements": 30,
            "n_species_with_usable_elements": 30,
            "n_selected_elements": 1088,
            "n_source_corrupt_exclusions": 2,
            "n_usable_elements": 1086,
            "exact_tesorter_join": True,
            "zero_depth_positions_retained": True,
            "expensive_upstream_tools_executed": False,
        }
        for key, expected in expected_manifest_values.items():
            self.assertEqual(manifest.get(key), expected, key)

        for record in manifest["data_outputs"].values():
            path = ROOT / record["path"]
            require_file(self, path)
            self.assertEqual(sha256(path), record["sha256"])
            self.assertEqual(len(csv_rows(path)), int(record["rows"]))
        for record in manifest["figures"]:
            path = ROOT / record["path"]
            require_file(self, path)
            self.assertEqual(sha256(path), record["sha256"])

        runner = SOURCE_R_RUNNER.read_text(encoding="utf-8")
        self.assertIn("ectopic_ltr30", runner)
        self.assertIn("ectopic_element_metrics_ltr30_v1.csv", runner)
        self.assertIn("ratio_terminal_internal_all_positions", runner)
        self.assertNotIn(HISTORICAL_LTR30_BASENAME, runner)

        review_manifest = json.loads(FIGURE_MANIFEST.read_text(encoding="utf-8"))
        self.assertEqual(review_manifest["scopes"]["ltr"]["n_species"], 30)
        self.assertEqual(review_manifest["scopes"]["ltr"]["n_elements"], 1086)
        review_manifest_text = json.dumps(review_manifest).lower()
        self.assertIn("zero", review_manifest_text)
        self.assertIn("all-position", review_manifest_text)
        self.assertIn("ectopic_element_metrics_ltr30_v1.csv", review_manifest_text)

        ltr_notebook = read_notebook(self, "04_ltr_deletion_footprint.ipynb")
        ltr_source = notebook_source(ltr_notebook)
        self.assertIn("ectopic_element_metrics_ltr30_v1.csv", ltr_source)
        self.assertIn("ectopic_excluded_elements_ltr30_v1.csv", ltr_source)
        self.assertNotIn(HISTORICAL_LTR30_BASENAME, ltr_source)
        frozen_sources = {row["source_path"] for row in csv_rows(FROZEN_REGISTRY)}
        self.assertFalse(
            any(Path(path).name == HISTORICAL_LTR30_BASENAME for path in frozen_sources)
        )

    def test_ltr_statistics_use_consistent_raw_log10_and_nonparametric_labels(self) -> None:
        stats_path = DATA_DIR / "ltr_terminal_internal_stats_v1.txt"
        require_file(self, stats_path)
        stats = stats_path.read_text(encoding="utf-8").lower()
        expected_labels = (
            r"raw(?:\s+ratio)?\s+one-way\s+anova",
            r"log10(?:[- ]transformed)?(?:\s+ratio)?\s+one-way\s+anova",
            r"kruskal-wallis(?:\s+test)?(?:\s+on)?\s+raw",
        )
        for pattern in expected_labels:
            self.assertRegex(stats, pattern)
        self.assertRegex(stats, r"unfiltered.{0,80}positive.{0,40}finite|positive.{0,40}finite.{0,80}unfiltered")
        self.assertRegex(stats, r"no\s+iqr\s+(?:filter|deletion|exclusion)")

        ltr_notebook = read_notebook(self, "04_ltr_deletion_footprint.ipynb")
        ltr_source = notebook_source(ltr_notebook).lower()
        for phrase in ("raw one-way anova", "log10 one-way anova", "kruskal-wallis raw"):
            self.assertIn(phrase, ltr_source)
        self.assertIn("ltr_terminal_internal_stats_v1.txt", ltr_source)

    def test_ltr_primary_view_screens_only_audited_mapping_artifact(self) -> None:
        require_file(self, LTR_ARTIFACT_SCREEN)
        rows = csv_rows(LTR_ARTIFACT_SCREEN)
        self.assertEqual(len(rows), 1086)
        self.assertEqual(len({row["species"] for row in rows}), 30)

        excluded = [row for row in rows if not normalized_bool(row["artifact_screen_keep"])]
        retained = [row for row in rows if normalized_bool(row["artifact_screen_keep"])]
        self.assertEqual(len(excluded), 1)
        self.assertEqual(len(retained), 1085)
        self.assertEqual(
            excluded[0]["element_id"],
            "JAUEJH010597481.1_De_5169_11831",
        )
        self.assertEqual(excluded[0]["artifact_screen_rule_id"], "one_sided_terminal_pileup_v1")
        self.assertGreater(abs(float(excluded[0]["left_right_terminal_log2_imbalance"])), 6.0)
        self.assertGreater(
            float(excluded[0]["ratio_terminal_internal_all_positions"]),
            float(excluded[0]["species_iqr_upper"]),
        )
        self.assertGreater(float(excluded[0]["fraction_of_species_ratio_sum"]), 0.5)
        self.assertLess(max(float(row["ratio_terminal_internal_all_positions"]) for row in retained), 14.0)
        self.assertEqual(len({row["species"] for row in retained}), 30)

        manifest = json.loads(FIGURE_MANIFEST.read_text(encoding="utf-8"))
        screen = manifest["ltr_mapping_artifact_screen"]
        self.assertEqual(screen["n_input_elements"], 1086)
        self.assertEqual(screen["n_flagged_elements"], 1)
        self.assertEqual(screen["n_retained_elements"], 1085)
        self.assertIs(screen["unfiltered_input_modified"], False)
        self.assertEqual(screen["left_right_abs_log2_threshold"], 6)

        notebook = read_notebook(self, "04_ltr_deletion_footprint.ipynb")
        source = notebook_source(notebook)
        self.assertIn("ltr_mapping_artifact_screen_ltr30_v1.csv", source)
        self.assertIn("ltr_terminal_internal_artifact_screened_ltr30_v1.png", source)
        self.assertRegex(source, r"(?i)unfiltered.{0,120}(?:retained|audit)")

    def test_selected_cluster_k_is_the_maximum_silhouette_not_hardcoded_five(self) -> None:
        metrics_path = DATA_DIR / "te_pca_cluster_metrics_te34_v1.csv"
        require_file(self, metrics_path)
        with metrics_path.open(newline="", encoding="utf-8-sig") as handle:
            rows = list(csv.DictReader(handle))
        self.assertGreaterEqual(len(rows), 2)
        self.assertTrue(
            {"k", "wss", "mean_silhouette", "selected"}.issubset(rows[0]),
            "cluster metrics columns",
        )
        selected = [row for row in rows if normalized_bool(row["selected"])]
        self.assertEqual(len(selected), 1, "exactly one cluster solution must be selected")
        selected_score = float(selected[0]["mean_silhouette"])
        maximum_score = max(float(row["mean_silhouette"]) for row in rows)
        self.assertAlmostEqual(selected_score, maximum_score, places=12)

        runner_source = SOURCE_R_RUNNER.read_text(encoding="utf-8")
        forbidden_hardcodes = (
            r"(?i)selected(?:_k)?\s*<-\s*5\b",
            r"(?i)selected(?:_k)?\s*=\s*5\b",
            r"(?i)kmeans\s*\([^\n]*centers\s*=\s*5\b",
        )
        for pattern in forbidden_hardcodes:
            self.assertIsNone(re.search(pattern, runner_source), pattern)

        repeat_notebook = read_notebook(self, "03_repeat_analysis_te34.ipynb")
        repeat_code = notebook_source(repeat_notebook, "code")
        self.assertIsNone(
            re.search(r"(?i)selected(?:_k)?\s*=\s*5\b", repeat_code),
            "notebook hardcodes k=5",
        )


if __name__ == "__main__":
    unittest.main()
