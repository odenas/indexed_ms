
MS_EXEC      = ../../../../fast_ms/bin/matching_stats_parallel.x
BS           = ../../block_sample.py
INPUT_FNAME  = ../Pan_troglodytes.CHIMP2.1.4.dna.chromosome.21.juststring
CHR          = 21

SAMPLE_SIZES = b32_n1024 b64_n512 b128_n256 b256_n128 b512_n64 b1024_n32 b32706069_n1
ANSWERS      = $(addsuffix .answer,$(SAMPLE_SIZES))
TARGETS      = $(addsuffix .hist,$(SAMPLE_SIZES))

all : $(TARGETS)


b32_n1024 : $(BS) $(INPUT_FNAME)
	python $^ 32 1024 >$@

b64_n512 : $(BS) $(INPUT_FNAME)
	python $^ 64 512 >$@

b128_n256 : $(BS) $(INPUT_FNAME)
	python $^ 128 256 >$@

b256_n128 : $(BS) $(INPUT_FNAME)
	python $^ 256 128 >$@

b512_n64 : $(BS) $(INPUT_FNAME)
	python $^ 512 64 >$@

b1024_n32 : $(BS) $(INPUT_FNAME)
	python $^ 1024 32 >$@

$(ANSWERS) : %.answer : % index 
	$(MS_EXEC) -s_path index -t_path $< -load_cst 1 -answer 1 >$@

$(TARGETS) : %.hist : %.answer
	cat $< \
		| sed -E -e 's/[[:blank:]]+/\n/g' \
		| awk '{A[$$1] += 1} END {for (i in A) print i, A[i]}' >$@

clean :
	rm -f $(ANSWERS) \
		$(subst answers,ms,$(ANSWERS)) $(subst answers,ms0,$(ANSWERS)) \
		$(subst answers,runs,$(ANSWERS)) $(subst answers,runs0,$(ANSWERS))
