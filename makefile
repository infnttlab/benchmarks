MAKE=make
MODULES=pi_loop primes
TARGETS=all clean
all:
	mkdir -p bin
	for b in $(MODULES) ; do $(MAKE) -C $$b $@ ; cp $$b/bin/* bin; done
clean:
	for b in $(MODULES) ; do $(MAKE) -C $$b $@ ; done
	rm -fr bin

