CC=	g++

BIN = $(DESTDIR)/usr/bin
FLAGS= -DMAXHITS=1000 -DTHREAD -funroll-loops -Lsamtools -Isamtools -Lgzstream -Igzstream -O3 -m64
#FLAGS= -DMAXHITS=1000 -funroll-loops -Lsamtools -Isamtools -Lgzstream -Igzstream -g -m64
#FLAGS= -DMAXHITS=1000 -funroll-loops -Lsamtools -Isamtools -Lgzstream -Igzstream -O3 -Wall -Wno-strict-aliasing -m64


THREAD=	-lpthread

SOURCE = align dbseq main pairs param reads utilities
OBJS1= $(patsubst %,%.o,$(SOURCE))

all: bsmap
%.o:%.cpp
	$(CC) $(FLAGS) -c $< -o $@
bsmap: $(OBJS1)
	(cd samtools; make)
	(cd gzstream; make)
	$(CC) $(FLAGS) $^ -o $@ $(THREAD) -lbam -lz -lgzstream
	rm -f *.o

clean:
	rm -f *.o *~ bsmap
	(cd samtools; make clean)
	(cd gzstream; make clean)
install:
	install -d $(BIN)
	install ./bsmap $(BIN)
	install ./sam2bam.sh $(BIN)
	install ./methratio.py $(BIN)
