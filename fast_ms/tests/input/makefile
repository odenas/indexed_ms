
INPUT_IDS   = rep_100000s_1000t_ab_1_5 rep_100000s_1000t_ab_2_10 rep_100000s_1000t_abcd_2_10 \
			  rep_100000s_1000t rnd_200_1024 rnd_200_128 rnd_200_16384 rnd_200_2048 rnd_200_256 \
			  rnd_200_32768 rnd_200_32 rnd_200_4096 rnd_200_512 rnd_200_64 rnd_20_10
DUMP_CST    = ../../bin/dump_cst.x
DUMP_MAXREP = ../../bin/dump_maxrep.x

CST         = $(addsuffix .s.rev.stree,$(INPUT_IDS))
MAXREP      = $(addsuffix .s.rev.maxrep,$(INPUT_IDS))

all : $(MAXREP)

%.rev.maxrep : % %.rev.stree | $(DUMP_MAXREP)
	$| -s_path $< -load_cst 1

$(CST) : %.rev.stree : % | $(DUMP_CST)
	$| -s_path $< 


clean :
	rm -f  $(CST) $(subst .rev.stree,.fwd.stree,$(CST)) \
		$(MAXREP) *.t.ms* *.t.runs*
