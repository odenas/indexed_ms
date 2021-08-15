

WRAPPERS = ms ms_slow range_queries range_queries_profile compress_ms extract_fasta \
		   diff_compress_ms diff_create_table \
		   longdiff_compress_ms longdiff_create_table \
		   nodiff_compress_ms nodiff_create_table

.PHONY : $(WRAPPERS) all
all : $(WRAPPERS)

$(WRAPPERS) :
	bash call_snakemake.sh $@

