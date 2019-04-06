
IID         = mm
DUMP_RIDX   = ../../fast_ms/bin/dump_range_index.x
COMPRESS_MS = ../../fast_ms/bin/compress_ms.x
RANGE_Q     = ../../fast_ms/bin/range_queries_profile.x -from_max_idx 80000000
RANGE_Q     = ../../fast_ms/bin/range_queries_profile.x -from_max_idx 800


RANGE_SIZES = $(addsuffix 0000000,2 4 6 10 20)
BLOCK_SIZES = 0 4 8 16 32 64 128 256 512 1024 2048 4096
BLOCK_SIZES = 0 4 128 512 1024 4096
BLOCK_SIZES = 0 32 1024

RANGES   = $(addsuffix .ridx,$(addprefix $(IID).t.ms.,$(BLOCK_SIZES)))
TARGETS  = $(addsuffix .csv,$(addprefix $(IID).t.ms.,$(BLOCK_SIZES)))

.SECONDARY : $(RANGES)

all :  $(TARGETS)

$(IID).t.ms.0.csv : $(IID).t.ms.none $(IID).t.ms.rrr $(IID).t.ms.0.ridx
	echo "compression,block_size,range_size,nqueries,method,time_ms" >$@
	for i in $(RANGE_SIZES); do \
		$(RANGE_Q) -ms_path $(IID).t.ms.none -compression none -ridx_path $< -range_size $$i -block_size 0 -niter 1 >> $@;\
		$(RANGE_Q) -ms_path $(IID).t.ms.rrr -compression rrr -ridx_path $< -range_size $$i -block_size 0 -niter 1 >> $@;\
	done

$(IID).t.ms.%.csv :  $(IID).t.ms.%.ridx
	echo "block_size,range_size,nqueries,method,time_ms" >$@
	for i in $(RANGE_SIZES); do \
		$(RANGE_Q) -ms_path $(IID).t.ms.none -compression none -ridx_path $< -range_size $$i -block_size $(subst $(IID).t.ms.,,$(subst .csv,,$@)) -niter 1000000 >> $@;\
		$(RANGE_Q) -ms_path $(IID).t.ms.rrr  -compression rrr  -ridx_path $< -range_size $$i -block_size $(subst $(IID).t.ms.,,$(subst .csv,,$@)) -niter 1000000 >> $@;\
	done

$(IID).t.ms.0.ridx : $(IID).t.ms
	touch $@

$(IID).t.ms.%.ridx : $(IID).t.ms
	$(DUMP_RIDX) -ms_path $<  -block_size $(subst $(IID).t.ms.,,$(subst .ridx,,$@)) -time_usage 1 >>$(IID).csv

$(IID).t.ms.rrr : $(IID).t.ms
	$(COMPRESS_MS) -ms_path $< -compression rrr

$(IID).t.ms.none : $(IID).t.ms
	cp $< $@

clean : 
	rm -f  $(TARGETS) $(IID).t.ms.csv
