
DIRS = box6 cluster6 feedback io kh orbit rt sedov sod zeldovich

.PHONY: tests clean

tests: $(DIRS)
	-for d in $(DIRS); do (cd $$d; $(MAKE) ); done

clean:
	-for d in $(DIRS); do (cd $$d; $(MAKE) clean ); done


