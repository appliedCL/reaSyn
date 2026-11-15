.PHONY: check-fibril-environment download-fibril-structures test analyse

PYTHON ?= .venv/bin/python

test:
	PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 $(PYTHON) -m pytest -q

analyse:
	$(PYTHON) -m reasyn.fibril_pipeline

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
