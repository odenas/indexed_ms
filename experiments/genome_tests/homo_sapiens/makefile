
DUMP_CST    = ../../../fast_ms/bin/dump_cst.x
SPECIES     = Homo_sapiens
CHROMOSOMES = 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y
CHR_FILES   = $(addsuffix .fa,$(addprefix $(SPECIES).GRCh38.dna.chromosome.,$(CHROMOSOMES)))
Z_CHR_FILES = $(addsuffix .gz,$(CHR_FILES))
STR_FILES   = $(subst .fa,.juststring,$(CHR_FILES))
CST_FILES   = $(addsuffix .fwd.stree,$(STR_FILES))

.SECONDARY : $(CHR_FILES) $(STR_FILES)

all : $(CST_FILES) $(SPECIES).juststring.fwd.stree


$(SPECIES).juststring.fwd.stree : $(SPECIES).juststring $(DUMP_CST)
	$(DUMP_CST) -s_path $< 

$(SPECIES).juststring : $(STR_FILES)
	cat $^ >$@

$(CST_FILES) : %.juststring.fwd.stree : %.juststring $(DUMP_CST)
	$(DUMP_CST) -s_path $< 

%.juststring : %.fa ../fa_to_juststring.py
	python ../fa_to_juststring.py --input $<  --output $@

#$(CHR_FILES) : %.fa : %.fa.gz
#	gzip -d $< 
#
#$(Z_CHR_FILES) : %.gz :
#	wget ftp://ftp.ensembl.org/pub/release-87/fasta/homo_sapiens/dna/$@


