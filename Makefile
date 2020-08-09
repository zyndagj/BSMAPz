CXX=	g++
SHELL=bash

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

.PHONY: conda
conda:
	cd conda && conda build --python 2.7 -c bioconda -c conda-forge -c defaults .

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
OS1=test_data/original_single.bsp.mr
OS2=test_data/original_single.sam.mr
OP=test_data/original_paired.sam.mr

.SILENT:
.PRECIOUS:
test_data/simulated.fastq.gz: test_data/simulated.fastq
	gzip -c $< > $@
test_data/test_%.sam.bam: test_data/test_%.sam
	samtools view -bS $< | samtools sort -o $@ 2> $@.log
	samtools index $@ 2>> $@.log
	@echo OK - converted $< to sorted and indexed BAM
# Test single end input
test_data/test_single.bsp test_data/test_single.sam test_data/test_single.bam: test_data/simulated.fastq | bsmapz
	./bsmapz -a $< -z 33 -p 2 -q 20 -d $(REF) -S 77345 -w 1000 -o $@ 2> $@.log
	@echo OK - Finished aligning $@
test_data/test_single%.mr: test_data/test_single%
	python methratio.py -z -r -d $(REF) -o $@ $< &> $@.log
	@echo OK - Finished calling methylation in $@
	diff -q $(OS1) $@
	@echo OK - $@ matches $(OS1)
	diff -q $(OS2) $@
	@echo OK - $@ matches $(OS2)
# Test compressed single end input
test_data/test_single_compressed.bsp test_data/test_single_compressed.sam test_data/test_single_compressed.bam: test_data/simulated.fastq.gz | bsmapz
	./bsmapz -a $< -z 33 -p 2 -q 20 -d $(REF) -S 77345 -w 1000 -o $@ 2> $@.log
	@echo OK - Finished aligning $@
# Test paired end input
test_data/test_paired.bsp test_data/test_paired.sam test_data/test_paired.bam: | test_data/simulated_1.fastq test_data/simulated_2.fastq bsmapz
	./bsmapz -a test_data/simulated_1.fastq -b test_data/simulated_2.fastq -z 33 -p 2 -q 20 -d $(REF) -S 77345 -w 1000 -o $@ 2> $@.log
	@echo OK - Finished aligning $@
test_data/test_paired.%.mr: test_data/test_paired.%
	python methratio.py -z -r -d $(REF) -o $@ $< &> $@.log || { cat $@.log; exit 1; }
	@echo OK - Finished calling methylation in $@
	diff -q $(OP) $@
	@echo OK - $@ matches $(OP)
test_data/test_paired_combine.%.mr: test_data/test_paired.%
	python methratio.py -g -z -r -d $(REF) -o $@ $< &> $@.log || { cat $@.log; exit 1; }
	@echo OK - Finished calling methylation in $@
	#diff -q $(OP) $@
	#@echo OK - $@ matches $(OP)
test_data/test_paired_type.%.mr: test_data/test_paired.%
	python methratio.py -x CG -z -r -d $(REF) -o $@ $< &> $@.log || { cat $@.log; exit 1; }
	@echo OK - Finished calling methylation in $@
	#diff -q $(OP) $@
	#@echo OK - $@ matches $(OP)
test_data/test_paired_ct.%.mr: test_data/test_paired.%
	python methratio.py -i no-action -z -r -d $(REF) -o $@ $< &> $@.log || { cat $@.log; exit 1; }
	@echo OK - Finished calling methylation in $@
	diff -q test_data/original_paired_ct.bam.mr $@
	@echo OK - $@ matches test_data/original_paired_ct.bam.mr
test_data/test_paired_int.bam.mr: test_data/test_paired.bam
	python methratio.py -z -I -r -d $(REF) -o $@ $< &> $@.log || { cat $@.log; exit 1; }
	@echo OK - Finished calling methylation in $@
	diff -q $(OP) $@
	@echo OK - $@ matches $(OP)
test_data/test_sam2bam.sam: test_data/test_paired.sam
	cd $(dir $<) && ln -s $(notdir $<) $(notdir $@)
test_data/test_sam2bam.bam: test_data/test_sam2bam.sam
	bash sam2bam.sh $< &> $@.log
	@echo OK - Finished testing sam2bam.sh

MR = $(shell echo test_data/test_{paired,single,single_compressed}.{sam.bam,bam,sam,bsp}.mr test_data/test_paired_{ct,type,combine,int}.bam.mr test_data/test_sam2bam.bam)
test: | bsmapz $(MR)
test-clean:
	rm -f test_data/test_{single,paired,sam2bam}*

clean:
	rm -f *.o *~ bsmapz
	(cd samtools; make clean)
	(cd gzstream; make clean)
install:
	install -d $(BIN)
	install ./bsmapz $(BIN)
	install ./sam2bam.sh $(BIN)
	install ./methratio.py $(BIN)
	install ./methdiff.py $(BIN)
