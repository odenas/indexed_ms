
NCPU = 4
RANGE_QUERIES = max_range_queries sum_range_queries


all: $(RANGE_QUERIES) split_ms

.PHONY: split_ms
split_ms: split_ms/Snakefile | ms_construction
	snakemake -s $< -j1 -d $(dir $<) --delete-all-output && snakemake -s $< -d $(dir $<) -j 4 -pr

.PHONY: $(RANGE_QUERIES)
$(RANGE_QUERIES) : %_range_queries: %_range_queries/Snakefile | ms_construction
	snakemake -s $< -j1 -d $(dir $<) --delete-all-output && snakemake -s $< -d $(dir $<) -j 4 -pr

.PHONY: ms_construction
ms_construction: ms_construction/Snakefile | input
	snakemake -s $< -j1 -d $(dir $<) --delete-all-output && snakemake -s $< -d $(dir $<) -j 1 -pr

.PHONY: input
input: input/Snakefile
	snakemake -s $< -d $(dir $<) -j $(NCPU)

clean:
	for d in input ms_construction $(RANGE_QUERIES); do \
		snakemake -s $$d/Snakefile -j1 -d $$d --delete-all-output; \
	done



