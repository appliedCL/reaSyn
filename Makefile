.PHONY: check-fibril-environment

check-fibril-environment:
	docker build --tag reasyn-fibril-environment .
	mkdir -p results
	docker run --rm --network none --mount 'type=bind,src=$(CURDIR)/results,dst=/workspace/results' reasyn-fibril-environment
	cp results/fibril_environment.json results/fibril_environment.first.json
	docker run --rm --network none --mount 'type=bind,src=$(CURDIR)/results,dst=/workspace/results' reasyn-fibril-environment
	cmp results/fibril_environment.first.json results/fibril_environment.json
	rm results/fibril_environment.first.json
