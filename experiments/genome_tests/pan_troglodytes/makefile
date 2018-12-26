
DUMP_CST    = ../../../fast_ms/bin/dump_cst.x
SPECIES     = Pan_troglodytes
CHROMOSOMES = 1 2A 2B 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y
Z_CHR_FILES = $(addsuffix .fa.gz,$(addprefix $(SPECIES).CHIMP2.1.4.dna.chromosome.,$(CHROMOSOMES)))
CHR_FILES   = $(addsuffix .fa,$(addprefix $(SPECIES).CHIMP2.1.4.dna.chromosome.,$(CHROMOSOMES)))
STR_FILES   = $(subst .fa,.juststring,$(CHR_FILES))
CST_FILES   = $(subst .juststring,.fwd.stree,$(STR_FILES))

.SECONDARY : $(CHR_FILES) $(STR_FILES)

all : $(CST_FILES)

%.fwd.stree : %.juststring $(DUMP_CST)
	$(DUMP_CST) -s_path $< 

%.juststring : %.fa ../fa_to_juststring.py
	python ../fa_to_juststring.py --input $<  --output $@

#$(CHR_FILES) : %.fa : %.fa.gz
#	gzip -d $< 
#
#$(Z_CHR_FILES) : %.gz :
#	wget ftp://ftp.ensembl.org/pub/release-87/fasta/homo_sapiens/dna/$@

