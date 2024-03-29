#define PRINT(name) \
if(d_flag[id##name]) { fprintf(f,"   %s\n",#name); d_flag[id##name] = 0; }
#define idAGN  0
#ifdef AGN
d_flag[idAGN] = 1;
#else
d_flag[idAGN] = 0;
#endif
#define idBLASTWAVE_FEEDBACK  1
#ifdef BLASTWAVE_FEEDBACK
d_flag[idBLASTWAVE_FEEDBACK] = 1;
#else
d_flag[idBLASTWAVE_FEEDBACK] = 0;
#endif
#define idCONSTANT_TIMESTEP  2
#ifdef CONSTANT_TIMESTEP
d_flag[idCONSTANT_TIMESTEP] = 1;
#else
d_flag[idCONSTANT_TIMESTEP] = 0;
#endif
#define idCOOLING  3
#ifdef COOLING
d_flag[idCOOLING] = 1;
#else
d_flag[idCOOLING] = 0;
#endif
#define idCOSMOLOGY  4
#ifdef COSMOLOGY
d_flag[idCOSMOLOGY] = 1;
#else
d_flag[idCOSMOLOGY] = 0;
#endif
#define idDEBUG_MEMORY_USE  5
#ifdef DEBUG_MEMORY_USE
d_flag[idDEBUG_MEMORY_USE] = 1;
#else
d_flag[idDEBUG_MEMORY_USE] = 0;
#endif
#define idDEBUG_MEMORY_USE_VERBOSE  6
#ifdef DEBUG_MEMORY_USE_VERBOSE
d_flag[idDEBUG_MEMORY_USE_VERBOSE] = 1;
#else
d_flag[idDEBUG_MEMORY_USE_VERBOSE] = 0;
#endif
#define idDENSITY_CHUNK_SIZE  7
#ifdef DENSITY_CHUNK_SIZE
d_flag[idDENSITY_CHUNK_SIZE] = 1;
#else
d_flag[idDENSITY_CHUNK_SIZE] = 0;
#endif
#define idDUST_EVOLUTION  8
#ifdef DUST_EVOLUTION
d_flag[idDUST_EVOLUTION] = 1;
#else
d_flag[idDUST_EVOLUTION] = 0;
#endif
#define idELECTRON_ION_NONEQUILIBRIUM  9
#ifdef ELECTRON_ION_NONEQUILIBRIUM
d_flag[idELECTRON_ION_NONEQUILIBRIUM] = 1;
#else
d_flag[idELECTRON_ION_NONEQUILIBRIUM] = 0;
#endif
#define idENRICHMENT  10
#ifdef ENRICHMENT
d_flag[idENRICHMENT] = 1;
#else
d_flag[idENRICHMENT] = 0;
#endif
#define idENRICHMENT_SNIa  11
#ifdef ENRICHMENT_SNIa
d_flag[idENRICHMENT_SNIa] = 1;
#else
d_flag[idENRICHMENT_SNIa] = 0;
#endif
#define idEXTRA_PRESSURE_SOURCE  12
#ifdef EXTRA_PRESSURE_SOURCE
d_flag[idEXTRA_PRESSURE_SOURCE] = 1;
#else
d_flag[idEXTRA_PRESSURE_SOURCE] = 0;
#endif
#define idGRAVITY  13
#ifdef GRAVITY
d_flag[idGRAVITY] = 1;
#else
d_flag[idGRAVITY] = 0;
#endif
#define idGRAVITY_IN_RIEMANN  14
#ifdef GRAVITY_IN_RIEMANN
d_flag[idGRAVITY_IN_RIEMANN] = 1;
#else
d_flag[idGRAVITY_IN_RIEMANN] = 0;
#endif
#define idHYDRO  15
#ifdef HYDRO
d_flag[idHYDRO] = 1;
#else
d_flag[idHYDRO] = 0;
#endif
#define idHYDRO_CHUNK_SIZE  16
#ifdef HYDRO_CHUNK_SIZE
d_flag[idHYDRO_CHUNK_SIZE] = 1;
#else
d_flag[idHYDRO_CHUNK_SIZE] = 0;
#endif
#define idHYDRO_TRACERS  17
#ifdef HYDRO_TRACERS
d_flag[idHYDRO_TRACERS] = 1;
#else
d_flag[idHYDRO_TRACERS] = 0;
#endif
#define idHYDRO_TRACERS_NGP  18
#ifdef HYDRO_TRACERS_NGP
d_flag[idHYDRO_TRACERS_NGP] = 1;
#else
d_flag[idHYDRO_TRACERS_NGP] = 0;
#endif
#define idINERT_GAS_TRACER  19
#ifdef INERT_GAS_TRACER
d_flag[idINERT_GAS_TRACER] = 1;
#else
d_flag[idINERT_GAS_TRACER] = 0;
#endif
#define idISOTROPIC_TURBULENCE_ENERGY  20
#ifdef ISOTROPIC_TURBULENCE_ENERGY
d_flag[idISOTROPIC_TURBULENCE_ENERGY] = 1;
#else
d_flag[idISOTROPIC_TURBULENCE_ENERGY] = 0;
#endif
#define idLOG_STAR_CREATION  21
#ifdef LOG_STAR_CREATION
d_flag[idLOG_STAR_CREATION] = 1;
#else
d_flag[idLOG_STAR_CREATION] = 0;
#endif
#define idMOMENTUM_DIFFUSION  22
#ifdef MOMENTUM_DIFFUSION
d_flag[idMOMENTUM_DIFFUSION] = 1;
#else
d_flag[idMOMENTUM_DIFFUSION] = 0;
#endif
#define idMPE_LOG  23
#ifdef MPE_LOG
d_flag[idMPE_LOG] = 1;
#else
d_flag[idMPE_LOG] = 0;
#endif
#define idMPI_MAX_MESSAGE_SIZE  24
#ifdef MPI_MAX_MESSAGE_SIZE
d_flag[idMPI_MAX_MESSAGE_SIZE] = 1;
#else
d_flag[idMPI_MAX_MESSAGE_SIZE] = 0;
#endif
#define idNDEBUG  25
#ifdef NDEBUG
d_flag[idNDEBUG] = 1;
#else
d_flag[idNDEBUG] = 0;
#endif
#define idnum_extra_energy_variables  26
#ifdef num_extra_energy_variables
d_flag[idnum_extra_energy_variables] = 1;
#else
d_flag[idnum_extra_energy_variables] = 0;
#endif
#define idOLDSTYLE_COOLING_EXPLICIT_SOLVER  27
#ifdef OLDSTYLE_COOLING_EXPLICIT_SOLVER
d_flag[idOLDSTYLE_COOLING_EXPLICIT_SOLVER] = 1;
#else
d_flag[idOLDSTYLE_COOLING_EXPLICIT_SOLVER] = 0;
#endif
#define idOLDSTYLE_PARTICLE_FILE_IGNORE_NGRID  28
#ifdef OLDSTYLE_PARTICLE_FILE_IGNORE_NGRID
d_flag[idOLDSTYLE_PARTICLE_FILE_IGNORE_NGRID] = 1;
#else
d_flag[idOLDSTYLE_PARTICLE_FILE_IGNORE_NGRID] = 0;
#endif
#define idOLDSTYLE_SF_ALGORITHM  29
#ifdef OLDSTYLE_SF_ALGORITHM
d_flag[idOLDSTYLE_SF_ALGORITHM] = 1;
#else
d_flag[idOLDSTYLE_SF_ALGORITHM] = 0;
#endif
#define idPAPI_PROFILING  30
#ifdef PAPI_PROFILING
d_flag[idPAPI_PROFILING] = 1;
#else
d_flag[idPAPI_PROFILING] = 0;
#endif
#define idPARTICLES  31
#ifdef PARTICLES
d_flag[idPARTICLES] = 1;
#else
d_flag[idPARTICLES] = 0;
#endif
#define idRADIATIVE_TRANSFER  32
#ifdef RADIATIVE_TRANSFER
d_flag[idRADIATIVE_TRANSFER] = 1;
#else
d_flag[idRADIATIVE_TRANSFER] = 0;
#endif
#define idREFINEMENT  33
#ifdef REFINEMENT
d_flag[idREFINEMENT] = 1;
#else
d_flag[idREFINEMENT] = 0;
#endif
#define idRT_ADD_EXTERNAL_LLS  34
#ifdef RT_ADD_EXTERNAL_LLS
d_flag[idRT_ADD_EXTERNAL_LLS] = 1;
#else
d_flag[idRT_ADD_EXTERNAL_LLS] = 0;
#endif
#define idRT_ADD_EXTERNAL_QSO_BACKGROUND  35
#ifdef RT_ADD_EXTERNAL_QSO_BACKGROUND
d_flag[idRT_ADD_EXTERNAL_QSO_BACKGROUND] = 1;
#else
d_flag[idRT_ADD_EXTERNAL_QSO_BACKGROUND] = 0;
#endif
#define idRT_BACKGROUND_HAARDT_MADAU  36
#ifdef RT_BACKGROUND_HAARDT_MADAU
d_flag[idRT_BACKGROUND_HAARDT_MADAU] = 1;
#else
d_flag[idRT_BACKGROUND_HAARDT_MADAU] = 0;
#endif
#define idRT_BACKGROUND_SELFCONSISTENT  37
#ifdef RT_BACKGROUND_SELFCONSISTENT
d_flag[idRT_BACKGROUND_SELFCONSISTENT] = 1;
#else
d_flag[idRT_BACKGROUND_SELFCONSISTENT] = 0;
#endif
#define idRT_CFI  38
#ifdef RT_CFI
d_flag[idRT_CFI] = 1;
#else
d_flag[idRT_CFI] = 0;
#endif
#define idRT_CHEMISTRY  39
#ifdef RT_CHEMISTRY
d_flag[idRT_CHEMISTRY] = 1;
#else
d_flag[idRT_CHEMISTRY] = 0;
#endif
#define idRT_CHEMISTRY_MINIMAL_MODEL  40
#ifdef RT_CHEMISTRY_MINIMAL_MODEL
d_flag[idRT_CHEMISTRY_MINIMAL_MODEL] = 1;
#else
d_flag[idRT_CHEMISTRY_MINIMAL_MODEL] = 0;
#endif
#define idRT_DUST_ABSORBS_BACKGROUND  41
#ifdef RT_DUST_ABSORBS_BACKGROUND
d_flag[idRT_DUST_ABSORBS_BACKGROUND] = 1;
#else
d_flag[idRT_DUST_ABSORBS_BACKGROUND] = 0;
#endif
#define idRT_DUST_CS  42
#ifdef RT_DUST_CS
d_flag[idRT_DUST_CS] = 1;
#else
d_flag[idRT_DUST_CS] = 0;
#endif
#define idRT_DUST_EVOLUTION  43
#ifdef RT_DUST_EVOLUTION
d_flag[idRT_DUST_EVOLUTION] = 1;
#else
d_flag[idRT_DUST_EVOLUTION] = 0;
#endif
#define idRT_EXACT_EOS  44
#ifdef RT_EXACT_EOS
d_flag[idRT_EXACT_EOS] = 1;
#else
d_flag[idRT_EXACT_EOS] = 0;
#endif
#define idRT_EXTERNAL_BACKGROUND  45
#ifdef RT_EXTERNAL_BACKGROUND
d_flag[idRT_EXTERNAL_BACKGROUND] = 1;
#else
d_flag[idRT_EXTERNAL_BACKGROUND] = 0;
#endif
#define idRT_FIXED_ISM  46
#ifdef RT_FIXED_ISM
d_flag[idRT_FIXED_ISM] = 1;
#else
d_flag[idRT_FIXED_ISM] = 0;
#endif
#define idRT_H2_RATE  47
#ifdef RT_H2_RATE
d_flag[idRT_H2_RATE] = 1;
#else
d_flag[idRT_H2_RATE] = 0;
#endif
#define idRT_HIGH_DENSITY  48
#ifdef RT_HIGH_DENSITY
d_flag[idRT_HIGH_DENSITY] = 1;
#else
d_flag[idRT_HIGH_DENSITY] = 0;
#endif
#define idRT_LWBANDS  49
#ifdef RT_LWBANDS
d_flag[idRT_LWBANDS] = 1;
#else
d_flag[idRT_LWBANDS] = 0;
#endif
#define idRT_LYMAN_ALPHA_HEATING  50
#ifdef RT_LYMAN_ALPHA_HEATING
d_flag[idRT_LYMAN_ALPHA_HEATING] = 1;
#else
d_flag[idRT_LYMAN_ALPHA_HEATING] = 0;
#endif
#define idRT_METHOD_OTVET  51
#ifdef RT_METHOD_OTVET
d_flag[idRT_METHOD_OTVET] = 1;
#else
d_flag[idRT_METHOD_OTVET] = 0;
#endif
#define idRT_NARROW_TABLE  52
#ifdef RT_NARROW_TABLE
d_flag[idRT_NARROW_TABLE] = 1;
#else
d_flag[idRT_NARROW_TABLE] = 0;
#endif
#define idRT_NO_TABLE  53
#ifdef RT_NO_TABLE
d_flag[idRT_NO_TABLE] = 1;
#else
d_flag[idRT_NO_TABLE] = 0;
#endif
#define idRT_OLDSTYLE_H2_SHIELDING  54
#ifdef RT_OLDSTYLE_H2_SHIELDING
d_flag[idRT_OLDSTYLE_H2_SHIELDING] = 1;
#else
d_flag[idRT_OLDSTYLE_H2_SHIELDING] = 0;
#endif
#define idRT_OTVET_CACHE_ET  55
#ifdef RT_OTVET_CACHE_ET
d_flag[idRT_OTVET_CACHE_ET] = 1;
#else
d_flag[idRT_OTVET_CACHE_ET] = 0;
#endif
#define idRT_OTVET_CACHE_RF  56
#ifdef RT_OTVET_CACHE_RF
d_flag[idRT_OTVET_CACHE_RF] = 1;
#else
d_flag[idRT_OTVET_CACHE_RF] = 0;
#endif
#define idRT_OTVET_NO_GLOBAL_ARRAY  57
#ifdef RT_OTVET_NO_GLOBAL_ARRAY
d_flag[idRT_OTVET_NO_GLOBAL_ARRAY] = 1;
#else
d_flag[idRT_OTVET_NO_GLOBAL_ARRAY] = 0;
#endif
#define idRT_OTVET_SAVE_FLUX  58
#ifdef RT_OTVET_SAVE_FLUX
d_flag[idRT_OTVET_SAVE_FLUX] = 1;
#else
d_flag[idRT_OTVET_SAVE_FLUX] = 0;
#endif
#define idRT_OUTPUT  59
#ifdef RT_OUTPUT
d_flag[idRT_OUTPUT] = 1;
#else
d_flag[idRT_OUTPUT] = 0;
#endif
#define idRT_PAH_CR  60
#ifdef RT_PAH_CR
d_flag[idRT_PAH_CR] = 1;
#else
d_flag[idRT_PAH_CR] = 0;
#endif
#define idRT_PARALLEL_USE_MPI  61
#ifdef RT_PARALLEL_USE_MPI
d_flag[idRT_PARALLEL_USE_MPI] = 1;
#else
d_flag[idRT_PARALLEL_USE_MPI] = 0;
#endif
#define idRT_TRACK_RECOMBINATION_PHOTONS  62
#ifdef RT_TRACK_RECOMBINATION_PHOTONS
d_flag[idRT_TRACK_RECOMBINATION_PHOTONS] = 1;
#else
d_flag[idRT_TRACK_RECOMBINATION_PHOTONS] = 0;
#endif
#define idRT_TRANSFER  63
#ifdef RT_TRANSFER
d_flag[idRT_TRANSFER] = 1;
#else
d_flag[idRT_TRANSFER] = 0;
#endif
#define idRT_TRANSFER_FLUX_CONSERVING  64
#ifdef RT_TRANSFER_FLUX_CONSERVING
d_flag[idRT_TRANSFER_FLUX_CONSERVING] = 1;
#else
d_flag[idRT_TRANSFER_FLUX_CONSERVING] = 0;
#endif
#define idRT_TRANSFER_METHOD  65
#ifdef RT_TRANSFER_METHOD
d_flag[idRT_TRANSFER_METHOD] = 1;
#else
d_flag[idRT_TRANSFER_METHOD] = 0;
#endif
#define idRT_UV  66
#ifdef RT_UV
d_flag[idRT_UV] = 1;
#else
d_flag[idRT_UV] = 0;
#endif
#define idRT_UV_OLDSTYLE_3x1  67
#ifdef RT_UV_OLDSTYLE_3x1
d_flag[idRT_UV_OLDSTYLE_3x1] = 1;
#else
d_flag[idRT_UV_OLDSTYLE_3x1] = 0;
#endif
#define idRT_XRAYS  68
#ifdef RT_XRAYS
d_flag[idRT_XRAYS] = 1;
#else
d_flag[idRT_XRAYS] = 0;
#endif
#define idSAVE_LOAD_BALANCE_PARTITION  69
#ifdef SAVE_LOAD_BALANCE_PARTITION
d_flag[idSAVE_LOAD_BALANCE_PARTITION] = 1;
#else
d_flag[idSAVE_LOAD_BALANCE_PARTITION] = 0;
#endif
#define idSF_FEEDBACK  70
#ifdef SF_FEEDBACK
d_flag[idSF_FEEDBACK] = 1;
#else
d_flag[idSF_FEEDBACK] = 0;
#endif
#define idSF_RECIPE  71
#ifdef SF_RECIPE
d_flag[idSF_RECIPE] = 1;
#else
d_flag[idSF_RECIPE] = 0;
#endif
#define idSTAR_FORMATION  72
#ifdef STAR_FORMATION
d_flag[idSTAR_FORMATION] = 1;
#else
d_flag[idSTAR_FORMATION] = 0;
#endif
#define idSTAR_PARTICLE_TYPES  73
#ifdef STAR_PARTICLE_TYPES
d_flag[idSTAR_PARTICLE_TYPES] = 1;
#else
d_flag[idSTAR_PARTICLE_TYPES] = 0;
#endif
#define idUNIQUE_RAND  74
#ifdef UNIQUE_RAND
d_flag[idUNIQUE_RAND] = 1;
#else
d_flag[idUNIQUE_RAND] = 0;
#endif
#define PRINT_ALL \
PRINT(AGN) \
PRINT(BLASTWAVE_FEEDBACK) \
PRINT(CONSTANT_TIMESTEP) \
PRINT(COOLING) \
PRINT(COSMOLOGY) \
PRINT(DEBUG_MEMORY_USE) \
PRINT(DEBUG_MEMORY_USE_VERBOSE) \
PRINT(DENSITY_CHUNK_SIZE) \
PRINT(DUST_EVOLUTION) \
PRINT(ELECTRON_ION_NONEQUILIBRIUM) \
PRINT(ENRICHMENT) \
PRINT(ENRICHMENT_SNIa) \
PRINT(EXTRA_PRESSURE_SOURCE) \
PRINT(GRAVITY) \
PRINT(GRAVITY_IN_RIEMANN) \
PRINT(HYDRO) \
PRINT(HYDRO_CHUNK_SIZE) \
PRINT(HYDRO_TRACERS) \
PRINT(HYDRO_TRACERS_NGP) \
PRINT(INERT_GAS_TRACER) \
PRINT(ISOTROPIC_TURBULENCE_ENERGY) \
PRINT(LOG_STAR_CREATION) \
PRINT(MOMENTUM_DIFFUSION) \
PRINT(MPE_LOG) \
PRINT(MPI_MAX_MESSAGE_SIZE) \
PRINT(NDEBUG) \
PRINT(num_extra_energy_variables) \
PRINT(OLDSTYLE_COOLING_EXPLICIT_SOLVER) \
PRINT(OLDSTYLE_PARTICLE_FILE_IGNORE_NGRID) \
PRINT(OLDSTYLE_SF_ALGORITHM) \
PRINT(PAPI_PROFILING) \
PRINT(PARTICLES) \
PRINT(RADIATIVE_TRANSFER) \
PRINT(REFINEMENT) \
PRINT(RT_ADD_EXTERNAL_LLS) \
PRINT(RT_ADD_EXTERNAL_QSO_BACKGROUND) \
PRINT(RT_BACKGROUND_HAARDT_MADAU) \
PRINT(RT_BACKGROUND_SELFCONSISTENT) \
PRINT(RT_CFI) \
PRINT(RT_CHEMISTRY) \
PRINT(RT_CHEMISTRY_MINIMAL_MODEL) \
PRINT(RT_DUST_ABSORBS_BACKGROUND) \
PRINT(RT_DUST_CS) \
PRINT(RT_DUST_EVOLUTION) \
PRINT(RT_EXACT_EOS) \
PRINT(RT_EXTERNAL_BACKGROUND) \
PRINT(RT_FIXED_ISM) \
PRINT(RT_H2_RATE) \
PRINT(RT_HIGH_DENSITY) \
PRINT(RT_LWBANDS) \
PRINT(RT_LYMAN_ALPHA_HEATING) \
PRINT(RT_METHOD_OTVET) \
PRINT(RT_NARROW_TABLE) \
PRINT(RT_NO_TABLE) \
PRINT(RT_OLDSTYLE_H2_SHIELDING) \
PRINT(RT_OTVET_CACHE_ET) \
PRINT(RT_OTVET_CACHE_RF) \
PRINT(RT_OTVET_NO_GLOBAL_ARRAY) \
PRINT(RT_OTVET_SAVE_FLUX) \
PRINT(RT_OUTPUT) \
PRINT(RT_PAH_CR) \
PRINT(RT_PARALLEL_USE_MPI) \
PRINT(RT_TRACK_RECOMBINATION_PHOTONS) \
PRINT(RT_TRANSFER) \
PRINT(RT_TRANSFER_FLUX_CONSERVING) \
PRINT(RT_TRANSFER_METHOD) \
PRINT(RT_UV) \
PRINT(RT_UV_OLDSTYLE_3x1) \
PRINT(RT_XRAYS) \
PRINT(SAVE_LOAD_BALANCE_PARTITION) \
PRINT(SF_FEEDBACK) \
PRINT(SF_RECIPE) \
PRINT(STAR_FORMATION) \
PRINT(STAR_PARTICLE_TYPES) \
PRINT(UNIQUE_RAND) \
;
