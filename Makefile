.PHONY: setup setup-python setup-r check-python check-env status catalog publication refresh identity validate test lint report figure-review figure-data figures methods-word te-bonsai te-pca path path-quick

setup: setup-python setup-r

setup-python:
	uv sync --locked --all-groups

setup-r:
	scripts/setup_r.sh

check-python:
	@test -x "$${UV_PROJECT_ENVIRONMENT:-.venv}/bin/python" || { echo 'Python environment missing; run make setup-python first.' >&2; exit 1; }
	uv run --locked --no-sync --offline --no-python-downloads python -B scripts/check_python_environment.py
	uv sync --locked --check --all-groups --offline --no-python-downloads

check-env: check-python
	scripts/setup_r.sh --check

status:
	uv run --locked --no-sync python -m scripts.project status

catalog:
	uv run --locked --no-sync python -m scripts.project catalog

publication:
	uv run --locked --no-sync python -m scripts.project publication
	uv run --locked --no-sync python -m scripts.project catalog

refresh: publication validate

identity:
	uv run --locked --no-sync python -m scripts.project identity

validate:
	uv run --locked --no-sync python -m scripts.project validate

test: validate
	uv run --locked --no-sync python -m pytest -q

lint:
	uv run --locked --no-sync ruff check scripts tests analyses --exclude '*.ipynb'

report: figures validate
	scripts/run_r.sh Rscript scripts/publication/render_study_report.R

figure-review: figures validate
	scripts/run_r.sh Rscript scripts/publication/render_figure_review.R

figure-data: validate
	uv run --locked --no-sync python scripts/publication/build_figure_data_workbook.py

methods-word:
	uv run --locked --no-sync python scripts/publication/build_methods_document.py --check

figures:
	scripts/run_r.sh Rscript scripts/publication/build_gbe_figures.R

te-bonsai:
	uv run --locked --no-sync python analyses/01_transposable_elements/bonsai_comparison/compare_bonsai_to_pca.py check
	uv run --locked --no-sync python analyses/01_transposable_elements/bonsai_comparison/validate_strategy_localization.py check

te-pca:
	uv run --locked --no-sync python analyses/01_transposable_elements/recompute_pca.py
	scripts/run_r.sh Rscript analyses/01_transposable_elements/analyze_pca_structure.R

path:
	uv run --locked --no-sync python analyses/03_phylogenetic_path/build_path_release.py
	scripts/run_r.sh Rscript analyses/03_phylogenetic_path/run_path_analysis.R

path-quick:
	uv run --locked --no-sync python analyses/03_phylogenetic_path/build_path_release.py
	scripts/run_r.sh Rscript analyses/03_phylogenetic_path/run_path_analysis.R --quick
