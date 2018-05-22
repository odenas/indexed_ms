
DUMP_CST        = /Users/denas/Library/Developer/Xcode/DerivedData/fast_ms-dtwaybjykudaehehgvtglnvhcjbp/Build/Products/Debug/dump_cst
BIN_TO_TXT      = /Users/denas/Library/Developer/Xcode/DerivedData/fast_ms-dtwaybjykudaehehgvtglnvhcjbp/Build/Products/Debug/bin_to_txt

MM = Mus_musculus
HS = Homo_sapiens
DR = Danio_rerio

MM_FA = $(wildcard $(MM)*.fa.gz)
HS_FA = $(wildcard $(HS)*.fa.gz)
DR_FA = $(wildcard $(DR)*.fa.gz)
ALL_FA = $(MM_FA) $(HS_FA) $(DR_FA)

.SECONDARY : $(addsuffix .juststring.fwd.stree,$(HS) $(MM) $(DR)) $(DR_FA)

all : $(addsuffix .hist,$(foreach P,$(MM_FA),$(subst .fa.gz,,$(P))__$(HS)) $(foreach P,$(HS_FA),$(subst .fa.gz,,$(P))__$(MM)) $(foreach P,$(HS_FA),$(subst .fa.gz,,$(P))__$(DR)))


# histograms

#/Users/denas/Library/Developer/Xcode/DerivedData/fast_ms-dtwaybjykudaehehgvtglnvhcjbp/Build/Products/Debug/bin_to_txt Homo_sapiens.GRCh38.dna.chromosome.10__Mus_musculus | awk '{A[$1] += 1} END {for(i in A) print i, A[i]}' 
%.hist : %
	$(BIN_TO_TXT) $< \
		| awk '{A[$$1] += 1} END {for(i in A) print i, A[i]}'  >$@

# MS vectors in binary
%__$(MM) : %.juststring $(MM).juststring.fwd.stree
	python fast_ms.py  --out_path $@ --load_cst --ms_progress 100 --runs_progress 100  --answer $(MM).juststring $< 

%__$(HS) : %.juststring $(HS).juststring.fwd.stree
	python fast_ms.py  --out_path $@ --load_cst --ms_progress 100 --runs_progress 100  --answer $(HS).juststring $< 

%__$(DR) : %.juststring $(DR).juststring.fwd.stree
	python fast_ms.py  --out_path $@ --load_cst --ms_progress 100 --runs_progress 100  --answer $(DR).juststring $< 

# suffix tree of a full species' genome 
%.fwd.stree : $(DUMP_CST) %
	$^

# raw sequence no Ns
$(MM).juststring : $(subst .fa.gz,.juststring,$(MM_FA))
	cat  $^ | tr -d '\n' >$@

$(HS).juststring : $(subst .fa.gz,.juststring,$(HS_FA))
	cat  $^ | tr -d '\n' >$@

$(DR).juststring : $(subst .fa.gz,.juststring,$(DR_FA))
	cat  $^ | tr -d '\n' >$@

%.juststring : %.fa.gz fa_to_juststring.py
	gzip -cd $< \
		| python fa_to_juststring.py --input /dev/stdin --output $@

# download sequences
Da%.fa.gz :  ; wget  ftp://ftp.ensembl.org/pub/release-87/fasta/danio_rerio/dna/$@

Ho%.fa.gz :  ; wget ftp://ftp.ensembl.org/pub/release-87/fasta/homo_sapiens/dna/$@

Mus%.fa.gz : ; wget ftp://ftp.ensembl.org/pub/release-87/fasta/mus_musculus/dna/$@

# some ucsc tools
faToTwoBit : 
	wget http://hgdownload.soe.ucsc.edu/admin/exe/macOSX.x86_64/faToTwoBit

twoBitToFa :
	wget http://hgdownload.soe.ucsc.edu/admin/exe/macOSX.x86_64/twoBitToFa
