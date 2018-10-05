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

##############################################
# Test section
##############################################
REF=test_data/test.fasta

test_data/simulated.fastq.gz: test_data/simulated.fastq
	gzip -c $< > $@
test_data/%.sam.bam: test_data/%.sam
	samtools view -bS $< | samtools sort -o $@ 2> $@.log
	samtools index $@ 2>> $@.log
	@echo OK - converted $< to sorted and indexed BAM
# Test single end input
test_data/single.%: test_data/simulated.fastq | bsmapz
	./bsmapz -a $< -z 33 -p 2 -q 20 -d $(REF) -S 77345 -w 1000 -o $@ 2> $@.log
	@echo OK - Finished aligning $@
test_data/single%.mr: test_data/single%
	python methratio.py -z -r -d $(REF) -o $@ $< 2> $@.log
	@echo OK - Finished calling methylation in $@
	diff -q original.mr $@
	@echo OK - $@ matches original.mr
# Test compressed single end input
test_data/single_compressed.%: test_data/simulated.fastq.gz | bsmapz
	./bsmapz -a $< -z 33 -p 2 -q 20 -d $(REF) -S 77345 -w 1000 -o $@ 2> $@.log
	@echo OK - Finished aligning $@
# Test paired end input
test_data/paired.%: | test_data/simulated_1.fastq test_data/simulated_2.fastq bsmapz
	./bsmapz -a test_data/simulated_1.fastq -b test_data/simulated_2.fastq -z 33 -p 2 -q 20 -d $(REF) -S 77345 -w 1000 -o $@ 2> $@.log
	@echo OK - Finished aligning $@
test_data/paired%.mr: test_data/paired%
	python methratio.py -z -r -d $(REF) -o $@ $< 2> $@.log
	@echo OK - Finished calling methylation in $@
	diff -q original_paired.mr $@
	@echo OK - $@ matches original_paired.mr

MR = $(shell test_data/{paired,single,single_compressed}.{sam.bam,bam,sam,bsp}.mr)
test: $(MR) | bsmapz

clean:
	rm -f *.o *~ bsmapz
	(cd samtools; make clean)
	(cd gzstream; make clean)
install:
	install -d $(BIN)
	install ./bsmapz $(BIN)
	install ./sam2bam.sh $(BIN)
	install ./methratio.py $(BIN)
