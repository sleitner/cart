
DIRS = cosmo_nbody  cosmo_adiabatic_hydro  cosmo_csf_hydro  cosmo_csf_hydro_rt  cosmo_csf_hydro_rt_full


TAGS = $(patsubst %,%/.phony,$(DIRS))

.PHONY: tests clean

tests: $(TAGS)

%/.phony:
	$(MAKE) --directory=$(*F)

clean:
	-for d in $(DIRS); do (cd $$d; $(MAKE) clean ); done


