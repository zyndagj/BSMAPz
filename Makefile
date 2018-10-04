CXX=	g++

# Install prefix
DESTDIR=/usr/local
# Install directory
BIN = $(DESTDIR)/bin
# Maximum number of equal best seeds to investigate
HITS=10000
# With pthreads
override FLAGS+= -DMAXHITS=$(HITS) -DTHREAD -funroll-loops -Lsamtools -Isamtools -Lgzstream -Igzstream -O3 -m64 -march=native
# Without pthreads
#override FLAGS+= -DMAXHITS=$(HITS) -funroll-loops -Lsamtools -Isamtools -Lgzstream -Igzstream -O3 -m64 -march=native

THREAD=	-lpthread

SOURCE = align dbseq main pairs param reads utilities
OBJS1= $(patsubst %,%.o,$(SOURCE))

all: bsmapz test

%.o:%.cpp
	$(CXX) $(FLAGS) -c $< -o $@

bsmapz: $(OBJS1)
	$(MAKE) -C samtools
	$(MAKE) -C gzstream
	$(CXX) $(FLAGS) $^ -o $@ $(THREAD) -lbam -lz -lgzstream

REF=test_data/test.fasta
test_data/simulated.fastq.gz: | test_data/simulated.fastq
	# Create compressed input
	gzip -c $| > $@
# Test aligner
test_data/simulated.sam: test_data/simulated.fastq bsmapz
	./bsmapz -a $< -z 33 -p 2 -q 20 -d $(REF) -S 77345 -w 10000 -o $@ 2> $@.log
	@echo Finished aligning $@
test_data/simulated.sam.bam: test_data/simulated.sam
	samtools view -bS $< | samtools sort -o $@
	samtools index $@
test_data/simulated.bam: test_data/simulated.fastq bsmapz
	./bsmapz -a $< -z 33 -p 2 -q 20 -d $(REF) -S 77345 -w 10000 -o $@ 2> $@.log
	@echo Finished aligning $@
test_data/simulated.bsp: test_data/simulated.fastq.gz bsmapz
	./bsmapz -a $< -z 33 -p 2 -q 20 -d $(REF) -S 77345 -w 10000 -o $@ 2> $@.log
	rm $<
	@echo Finished aligning $@
# Test methratio
test_data/simulated.%.mr: test_data/simulated.%
	python methratio.py -z -r -d $(REF) -o $@ $<
	rm $<
	@echo Finished creating $@


REF_MR = test_data/original.mr
MR = test_data/simulated.sam.mr test_data/simulated.bam.mr test_data/simulated.bsp.mr test_data/simulated.sam.bam.mr
.SILENT: test
test: bsmapz $(MR)
	@echo Comparing methratio outputs
	for a in $(MR) $(REF_MR); do \
		for b in $(MR) $(REF_MR); do \
			if [ $$a != $$b ]; then \
				diff -q $$a $$b && echo "OK   - $$a and $$b MATCH" || echo "FAIL - $$a and $$b DIFFER"; \
			fi; \
		done; \
	done; \

clean:
	rm -f *.o *~ bsmapz
	(cd samtools; make clean)
	(cd gzstream; make clean)
install:
	install -d $(BIN)
	install ./bsmapz $(BIN)
	install ./sam2bam.sh $(BIN)
	install ./methratio.py $(BIN)
