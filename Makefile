HDRS:=$(wildcard src/*.hh)
CTRS:=$(wildcard src/*.cc)
LIBS:=-lz

main: hl-csa-raptor.o
	@echo "Try it with 'make test'"

test: _data/London hl-csa-raptor.o
	./hl-csa-raptor.o -nq=5 $<

%.o: src/%.cc $(HDRS)
	g++ -std=c++11 -O3 -pthread -o $@ $< $(LIBS)


REPO:=https://files.inria.fr/gang/graphs/public_transport

_data/%:
	mkdir -p _data/$*
	cd _data/$*; \
	for f in queries-unif.csv stop_times.csv.gz transfers.csv.gz walk_and_transfer_inhubs.gr.gz walk_and_transfer_outhubs.gr.gz; do \
		curl -o $$f $(REPO)/$*/$$f ;\
	done

clean:
	rm -fr *~ */*~ _* *.o
