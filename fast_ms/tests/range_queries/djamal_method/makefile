
EXEC  = ../../../bin/range_queries.x
INPUT = ../../input/rnd_200_256.t.ms
END = $(shell python -c "print(' '.join([str(i) for i in range(1, 1000)]))")

all : baseline.txt djamal.txt
	diff $^


baseline.txt :
	for i in $(END); do \
		$(EXEC) -ms_path $(INPUT) -from_idx 0 -to_idx $$i -block_size 0 -compression none \
			| awk -F ": " '{print $$2}' >>$@;\
	done

djamal.txt :
	for i in $(END); do \
		$(EXEC) -ms_path $(INPUT) -from_idx 0 -to_idx $$i -block_size -1 -compression none \
			| awk -F ": " '{print $$2}' >>$@;\
	done

clean :
	rm -f baseline.txt djamal.txt
