
FAST_MS = /Users/denas/Library/Developer/Xcode/DerivedData/fast_ms-dtwaybjykudaehehgvtglnvhcjbp/Build/Products/Debug/fd_ms
SLOW_MS = slow_ms.py

.PHONY: test clean rtime



test: | script_tests.py $(FAST_MS) $(SLOW_MS) a.txt
	python $< $| 

rtime : | $(FAST_MS) $(SLOW_MS) a.txt
	time head -1  a.txt | $(FAST_MS) /dev/stdin
	time head -1  a.txt | python $(SLOW_MS) /dev/stdin

a.txt : generate_strings.py
	python $< --l 6000 10 >$@

clean  :
	rm -f a.txt
