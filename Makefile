.PHONY: setup setup-python setup-r status catalog publication refresh identity validate test lint report figure-review figure-data figures te-pca path path-quick

setup: setup-python setup-r

setup-python:
	uv sync --all-groups

setup-r:
	scripts/setup_r.sh

status:
	uv run python -m scripts.project status

catalog:
	uv run python -m scripts.project catalog

publication:
	uv run python -m scripts.project publication
	uv run python -m scripts.project catalog

refresh: publication validate

identity:
	uv run python -m scripts.project identity

validate:
	uv run python -m scripts.project validate

test: validate
	uv run pytest -q

lint:
	uv run ruff check scripts tests analyses --exclude '*.ipynb'

report: figures validate
	scripts/run_r.sh Rscript scripts/publication/render_study_report.R

figure-review: figures validate
	scripts/run_r.sh Rscript scripts/publication/render_figure_review.R

figure-data: validate
	uv run --frozen --no-sync python scripts/publication/build_figure_data_workbook.py

figures:
	scripts/run_r.sh Rscript scripts/publication/build_gbe_figures.R

te-pca:
	uv run python analyses/01_transposable_elements/recompute_pca.py
	scripts/run_r.sh Rscript analyses/01_transposable_elements/analyze_pca_structure.R

path:
	scripts/run_r.sh Rscript analyses/03_phylogenetic_path/run_path_analysis.R

path-quick:
	scripts/run_r.sh Rscript analyses/03_phylogenetic_path/run_path_analysis.R --quick
