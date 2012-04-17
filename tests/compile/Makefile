
DIRS = cosmo_nbody  cosmo_adiabatic_hydro  cosmo_csf_hydro  cosmo_csf_hydro_rt  cosmo_csf_hydro_rt_full plain_hydro plain_hydro_gravity


TAGS = $(patsubst %,%/.phony,$(DIRS))

.PHONY: tests clean

tests: $(TAGS)

%/.phony:
	$(MAKE) --directory=$(*F)

clean:
	-for d in $(DIRS); do (cd $$d; $(MAKE) clean ); done


