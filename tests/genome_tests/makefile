
SPECIES   = Mus_musculus Homo_sapiens Danio_rerio
HUMAN_CHR = $(wildcard genome_data/Ho*2.juststring)
MOUSE_CHR = $(wildcard genome_data/Mu*.juststring)
ZFISH_CHR = $(wildcard genome_data/Da*.juststring)


#Homo_sapiens.GRCh38.dna.chromosome.MT__Mus_musculus.mk.pkl

HUMAN_MOUSE_TARGETS = $(foreach C,$(HUMAN_CHR),$(addsuffix __Mus_musculus.mk.pkl,$(subst .juststring,,$(C))))

all : $(HUMAN_MOUSE_TARGETS)

Homo_sapiens.GRCh38.dna.chromosome.%__Mus_musculus.mk.pkl : 
	./script_benchmark.py $(subst __Mus_musculus.mk.pkl,,$@).juststring Mus_musculus --load_cst --runs_progress 100 --ms_progress 100


