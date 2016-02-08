CXX=	g++

BIN = $(DESTDIR)/usr/bin
# With pthreads
override FLAGS+= -DMAXHITS=1000 -DTHREAD -funroll-loops -Lsamtools -Isamtools -Lgzstream -Igzstream -O3 -m64 -march=native
# Without pthreads
#override FLAGS+= -DMAXHITS=1000 -funroll-loops -Lsamtools -Isamtools -Lgzstream -Igzstream -O3 -m64 -march=native

THREAD=	-lpthread

SOURCE = align dbseq main pairs param reads utilities
OBJS1= $(patsubst %,%.o,$(SOURCE))

all: bsmap
%.o:%.cpp
	$(CXX) $(FLAGS) -c $< -o $@
bsmap: $(OBJS1)
	(cd samtools; make)
	(cd gzstream; make)
	$(CXX) $(FLAGS) $^ -o $@ $(THREAD) -lbam -lz -lgzstream
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
