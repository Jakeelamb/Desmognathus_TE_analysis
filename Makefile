PIXI ?= pixi
PIXI_RUN = $(PIXI) run --manifest-path pixi.toml --environment default --clean-env --locked --no-install
LOCAL_TASKS = status catalog publication refresh identity validate test lint report figure-review figure-data figures methods-word te-bonsai te-pca path path-quick

.PHONY: setup setup-python setup-r check-python check-env _runtime $(LOCAL_TASKS)

setup:
	$(PIXI) install --manifest-path pixi.toml --locked --all
	$(PIXI_RUN) setup

setup-python setup-r: setup

_runtime:
	@test -x .pixi/envs/default/bin/python && test -x .pixi/envs/default/bin/Rscript || { echo 'Pixi environment missing; run make setup first.' >&2; exit 1; }

check-python check-env: _runtime
	$(PIXI_RUN) --offline $@

$(LOCAL_TASKS): _runtime
	$(PIXI_RUN) $@
