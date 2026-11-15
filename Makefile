.PHONY: check-fibril-environment download-fibril-structures test analyse notebook causal-overlap

PYTHON ?= .venv/bin/python

test:
	PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 $(PYTHON) -m pytest -q

analyse:
	$(PYTHON) -m reasyn.fibril_pipeline

notebook:
	$(PYTHON) -m jupyter execute notebooks/analyse_fibril_structures.ipynb --timeout=1800 --inplace

causal-overlap:
	docker run --rm --network none --mount 'type=bind,src=$(CURDIR)/data,dst=/workspace/data' --mount 'type=bind,src=$(CURDIR)/data/raw,dst=/workspace/data/raw,readonly' --mount 'type=bind,src=$(CURDIR)/notebooks,dst=/workspace/notebooks' reasyn-analysis jupyter execute notebooks/assess_fibril_causal_overlap.ipynb --timeout=120 --inplace

check-fibril-environment:
	docker build --tag reasyn-fibril-environment .
	mkdir -p results
	docker run --rm --network none --mount 'type=bind,src=$(CURDIR)/results,dst=/workspace/results' reasyn-fibril-environment
	cp results/fibril_environment.json results/fibril_environment.first.json
	docker run --rm --network none --mount 'type=bind,src=$(CURDIR)/results,dst=/workspace/results' reasyn-fibril-environment
	cmp results/fibril_environment.first.json results/fibril_environment.json
	rm results/fibril_environment.first.json

download-fibril-structures:
	docker build --tag reasyn-fibril-environment .
	mkdir -p data
	docker run --rm --mount 'type=bind,src=$(CURDIR)/data,dst=/workspace/data' --mount 'type=bind,src=$(CURDIR)/notebooks,dst=/workspace/notebooks' reasyn-fibril-environment jupyter execute notebooks/download_fibril_structures.ipynb --timeout=1800 --inplace
