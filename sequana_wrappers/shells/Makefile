MAKEFILES := $(shell find . -mindepth 3 -name Makefile)

clean:
	@for mk in $(MAKEFILES); do \
	    $(MAKE) -C $$(dirname $$mk) clean; \
	done
