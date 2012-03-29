
DIRS = box6 cluster6 io kh orbit rt sedov sod zeldovich

TAGS = $(patsubst %,%/.phony,$(DIRS))

.PHONY: tests clean

tests: $(TAGS)

%/.phony:
	$(MAKE) --directory=$(*F)

clean:
	-for d in $(DIRS); do (cd $$d; $(MAKE) clean ); done


