.PHONY: setup setup-python setup-r status catalog publication refresh identity validate test lint notebooks figures te-pca te-pca-view path path-quick

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

notebooks:
	uv run jupyter lab analyses

figures:
	scripts/run_r.sh Rscript scripts/publication/build_gbe_figures.R

te-pca:
	uv run python analyses/01_transposable_elements/recompute_pca.py
	scripts/run_r.sh Rscript analyses/01_transposable_elements/analyze_pca_structure.R

te-pca-view:
	uv run jupyter lab analyses/01_transposable_elements/explore.ipynb

path:
	scripts/run_r.sh Rscript analyses/03_phylogenetic_path/run_path_analysis.R

path-quick:
	scripts/run_r.sh Rscript analyses/03_phylogenetic_path/run_path_analysis.R --quick
