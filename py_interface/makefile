
all : install

.PHONY: install
install:
	pip install --verbose -e .
	ln -s -f ./src mstat

.PHONY: clean
clean :
	pip uninstall mstat
