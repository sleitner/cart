
DIRS = et  iliev0  iliev1  iliev1a  iliev5  iliev6

TAGS = $(patsubst %,%/.phony,$(DIRS))

.PHONY: tests clean

tests: $(TAGS)

%/.phony:
	$(MAKE) --directory=$(*F)

clean:
	-for d in $(DIRS); do (cd $$d; $(MAKE) clean ); done


