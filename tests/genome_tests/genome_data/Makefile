
DUMP_CST        = /Users/denas/Library/Developer/Xcode/DerivedData/fast_ms-dtwaybjykudaehehgvtglnvhcjbp/Build/Products/Debug/dump_cst

MM_FTP_PATH     = ftp://ftp.ensembl.org/pub/release-87/fasta/mus_musculus/dna
HS_FTP_PATH     = ftp://ftp.ensembl.org/pub/release-87/fasta/homo_sapiens/dna
DR_FTP_PATH     = ftp://ftp.ensembl.org/pub/release-87/fasta/danio_rerio/dna

MM_CHROMOSOMES  = 1 10 11 12 13 14 15 16 17 18 19 2 3 4 5 6 7 8 9 MT X Y 10
HS_CHROMOSOMES  = 1 10 11 12 13 14 15 16 17 18 19 2 20 21 22 3 4 5 6 7 8 9 MT X Y
DR_CHROMOSOMES  = 1 10 11 12 13 14 15 16 17 18 19 2 20 21 22 23 24 25 3 4 5 6 7 8 9 MT

MM_FASTA        = $(addsuffix .fa,$(addprefix Mus_musculus.GRCm38.dna.chromosome.,$(MM_CHROMOSOMES)))
HS_FASTA        = $(addsuffix .fa,$(addprefix Homo_sapiens.GRCh38.dna.chromosome.,$(HS_CHROMOSOMES)))
DR_FASTA        = $(addsuffix .fa,$(addprefix Danio_rerio.GRCz10.dna.chromosome.,$(DR_CHROMOSOMES)))


#all : homo_sapiens.GRCh38.87.2bit mus_musculus.GRCm38.87.2bit danio_rerio.GRCz10.87.2bit

all : $(subst .fa,.fwd.stree,$(MM_FASTA)) $(subst .fa,.fwd.stree,$(HS_FASTA)) $(subst .fa,.fwd.stree,$(DR_FASTA))

Da%.fa.gz :  ; wget  ftp://ftp.ensembl.org/pub/release-87/fasta/danio_rerio/dna/$@

Ho%.fa.gz :  ; wget ftp://ftp.ensembl.org/pub/release-87/fasta/homo_sapiens/dna/$@

Mus%.fa.gz : ; wget ftp://ftp.ensembl.org/pub/release-87/fasta/mus_musculus/dna/$@

homo_sapiens.GRCh38.87.2bit : $(HS_FASTA) | faToTwoBit
	./faToTwoBit $^ $@

mus_musculus.GRCm38.87.2bit : $(MM_FASTA) | faToTwoBit
	./faToTwoBit $^ $@

danio_rerio.GRCz10.87.2bit : $(DR_FASTA) | faToTwoBit
	./faToTwoBit $^ $@

%.fwd.stree : %.juststring
	$(DUMP_CST) $< 

%.juststring : %.fa
	python fa_to_juststring.py --input $< --output $@

%.fa : %.fa.gz
	gzip -d $<

faToTwoBit : 
	wget http://hgdownload.soe.ucsc.edu/admin/exe/macOSX.x86_64/faToTwoBit
