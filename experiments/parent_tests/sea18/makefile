
SUBDIRS = $(wildcard trial_*)

all :
	for i in  $(SUBDIRS) ; \
	do \
		$(MAKE) -C $$i; \
	done;

clean:
	for i in  $(SUBDIRS) ; \
	do \
		$(MAKE) -C $$i clean ; \
	done;

clean_data:
	for i in  $(SUBDIRS) ; \
	do \
		$(MAKE) -C $$i clean_data ; \
	done;

