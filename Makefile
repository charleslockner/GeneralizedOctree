CC=g++
EXE=test
CFLAGS=-std=c++11 -I. -O3 -g -DMACOSX -MMD

.PHONY: all run clean

all: $(EXE)

run: $(EXE)
	./$(EXE)

clean:
	rm -rf $(EXE) *.DS_Store *~

# Special rule for model.test (needs geometry)
$(EXE): test.cpp geometry.cpp octree.cpp
	$(CC) $(CFLAGS) -o $@ $^
