
FAST_MS = /Users/denas/Library/Developer/Xcode/DerivedData/fast_ms-dtwaybjykudaehehgvtglnvhcjbp/Build/Products/Debug/fd_ms
SLOW_MS = slow_ms.py

.PHONY: test



test: | script_tests.py $(FAST_MS) $(SLOW_MS) a.txt
	python $< $| 

a.txt : generate_strings.py
	python $< 30 >$@
