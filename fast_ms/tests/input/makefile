
__IMPORT_STR = import sys; sys.path.append('..'); from _paths import *
DUMP_CST     = $(shell python -c "$(__IMPORT_STR); print(dump_cst)")
DUMP_MAXREP  = $(shell python -c "$(__IMPORT_STR); print(dump_maxrep)")
MS_SLOW      = $(shell python -c "$(__IMPORT_STR); print(ms_slow)")
INPUT_IDS    = $(shell python -c "$(__IMPORT_STR); print(' '.join(iids))")


CST          = $(shell python -c "$(__IMPORT_STR); print(' '.join(ipair(i).rev_cst for i in iids))")
MAXREP       = $(shell python -c "$(__IMPORT_STR); print(' '.join(ipair(i).maxrep for i in iids))")
MSTAT        = $(shell python -c "$(__IMPORT_STR); print(' '.join(ipair(i).mstat for i in iids))")


all : $(MAXREP) $(MSTAT)

$(MAXREP) : %.rev.maxrep : % %.rev.stree | $(DUMP_MAXREP)
	$| -s_path $< -load_cst 1

$(CST) : %.rev.stree : % | $(DUMP_CST)
	$| -s_path $< 

$(MSTAT) : %.mstat : %.s %.t | $(MS_SLOW)
	$| -s_path $< -t_path $(subst .s,.t,$<) >$@

clean :
	rm -f  $(CST) $(subst .rev.stree,.fwd.stree,$(CST)) \
		$(MAXREP) *.t.ms* *.t.runs*
