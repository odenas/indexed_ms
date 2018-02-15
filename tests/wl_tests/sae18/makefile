
SUBDIRS = $(wildcard trial_*)

stats : data
	for i in  $(SUBDIRS) ; \
	do \
		$(MAKE) -C $$i -f makefile.stats; \
	done;

clean_stats :
	for i in  $(SUBDIRS) ; \
	do \
		$(MAKE) -C $$i -f makefile.stats clean ; \
	done;

data :
	for i in  $(SUBDIRS) ; \
	do \
		$(MAKE) -C $$i; \
	done;

clean_data :
	for i in  $(SUBDIRS) ; \
	do \
		$(MAKE) -C $$i clean ; \
	done;

