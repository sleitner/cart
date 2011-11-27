
DIRS = et  iliev0  iliev1  iliev1a  iliev1b  iliev5  iliev6

.PHONY: tests clean

tests: $(DIRS)
	-for d in $(DIRS); do (cd $$d; $(MAKE) ); done

clean:
	-for d in $(DIRS); do (cd $$d; $(MAKE) clean ); done


