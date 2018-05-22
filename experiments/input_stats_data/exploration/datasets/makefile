
BIN         = ../../../../fast_ms/bin
DUMP_CST    = $(BIN)/dump_cst.x
DUMP_MAXREP = $(BIN)/dump_maxrep.x
ISTATS      = ./input_stats.py

ODIR        = ./data_dir
I_SCRIPT    = ../../../datasets/generate_input3.py
S_SCRIPT    = gen_script.py
CFG         = input_spec.json

S_OUTPUTS   = $(shell python $(S_SCRIPT) $(CFG) $(ODIR) | awk '$$1=="mv" && NR % 2 == 0 {print $$4}')
ST_OUTPUTS  = $(addsuffix .fwd.stree,$(S_OUTPUTS))
MR_OUTPUTS  = $(addsuffix .rev.maxrep,$(S_OUTPUTS))
CSV_OUTPUTS = $(addsuffix .csv,$(S_OUTPUTS))

.SECONDARY : $(ST_OUTPUTS) $(MR_OUTPUTS)

all : $(CSV_OUTPUTS)


script.sh : $(S_SCRIPT) $(CFG) $(I_SCRIPT) | $(ODIR)
	python $< $(CFG) $(ODIR) --script $(I_SCRIPT) >$@

$(ODIR)/done : script.sh $(I_SCRIPT)
	bash $<

%.fwd.stree : $(DUMP_CST) $(ODIR)/done
	$(DUMP_CST) -s_path $(subst .fwd.stree,,$@)

%.rev.maxrep : %.fwd.stree $(DUMP_MAXREP)
	$(DUMP_MAXREP) -load_cst 1 -s_path $(subst .rev.maxrep,,$@)

%.csv  : $(ISTATS) $(ODIR)/done %.fwd.stree %.rev.maxrep
	$(ISTATS) --repeat 25 $(subst .csv,,$@) $(subst .s.csv,.t,$@) --load_maxrep --load_cst >$@

$(ODIR) : 
	mkdir $@

clean : 
	rm -rf $(ODIR) script.sh
