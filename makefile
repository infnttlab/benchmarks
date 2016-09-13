MAKE=make
MODULES=pi_loop primes
TARGETS=all clean

all:
	for b in $(MODULES) ; do $(MAKE) -C $$b $@ ; cp $$b/bin/* bin; done

clean:
	for b in $(MODULES) ; do $(MAKE) -C $$b $@; done
	rm -fr bin

