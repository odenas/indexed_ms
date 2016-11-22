
FAST_MS = /Users/denas/Library/Developer/Xcode/DerivedData/fast_ms-dtwaybjykudaehehgvtglnvhcjbp/Build/Products/Debug/fd_ms
SLOW_MS = slow_ms.py
BDIR    = input_data
PROTEIN_DSETS = $(addprefix proteins.50MB1000_,1000 5000 10000 20000 40000 80000 160000 320000)
N = 10

.PHONY: test clean input

protein_benchmark.csv :  script_benchmark.py $(FAST_MS)
	python $< --output $@ $(BDIR) $(PROTEIN_DSETS)

test:  script_tests.py $(SLOW_MS) $(FAST_MS)  | $(BDIR)
	python $^ $(N) --base_dir $| 

random_input : generate_input.py | $(BDIR)
	python $< --l 100  --L 120 --base_dir $| $(N)

clean  :
	rm -f a.txt






