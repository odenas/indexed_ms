
SPECIES   = Mus_musculus Homo_sapiens Danio_rerio
HUMAN_CHR = $(wildcard genome_data/Ho*.juststring)
MOUSE_CHR = $(wildcard genome_data/Mu*.juststring)
ZFISH_CHR = $(wildcard genome_data/Da*.juststring)


#Homo_sapiens.GRCh38.dna.chromosome.MT__Mus_musculus.mk.pkl

HUMAN_MOUSE_TARGETS = $(foreach C,$(HUMAN_CHR),$(addsuffix __Mus_musculus,$(subst .juststring,,$(C))))

all : $(HUMAN_MOUSE_TARGETS)
	echo $(HUMAN_CHR)

#denas@denas-osx:$ ./script_benchmark.py genome_data/Homo_sapiens.GRCh38.dna.chromosome.MT.juststring Mus_musculus --load_cst --runs_progress 5 --ms_progress 5

Homo_sapiens.GRCh38.dna.chromosome.%__Mus_musculus : 
	./script_benchmark.py $(subst __Mus_musculus,,$@).juststring Mus_musculus --load_cst --runs_progress 100 --ms_progress 100

Homo_sapiens.GRCh38.dna.chromosome.%__Danio_rerio :
	./script_benchmark.py $(subst __Danio_rerio,,$@).juststring Danio_rerio --load_cst --runs_progress 100 --ms_progress 100


