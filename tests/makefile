
FAST_MS = /Users/denas/Library/Developer/Xcode/DerivedData/fast_ms-dtwaybjykudaehehgvtglnvhcjbp/Build/Products/Debug/fd_ms
SLOW_MS = slow_ms.py
BDIR    = input_data
N = 10

.PHONY: test clean input

all : input test


test:  script_tests.py $(SLOW_MS) $(FAST_MS)  | $(BDIR)
	python $^ $(N) --base_dir $| 

input : generate_strings.py | $(BDIR)
	python $< --l 100  --L 120 --base_dir $| $(N)

clean  :
	rm -f a.txt
