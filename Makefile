TOP_DIR=.
include Make.config

.phony: all compile_tests run_tests clean docs utils

all:	compile_tests 

docs:
	cd docs; $(MAKE)

utils:
	cd utils; $(MAKE)

compile_tests:
	cd tests; $(MAKE)

run_tests:
	cd tests; $(MAKE) run

clean:
	cd tests; $(MAKE) clean
	cd utils; $(MAKE) clean
	cd docs; $(MAKE) clean
