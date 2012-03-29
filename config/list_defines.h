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
#define idDEBUG  5
#ifdef DEBUG
d_flag[idDEBUG] = 1;
#else
d_flag[idDEBUG] = 0;
#endif
#define idDEBUG_MEMORY_USE  6
#ifdef DEBUG_MEMORY_USE
d_flag[idDEBUG_MEMORY_USE] = 1;
#else
d_flag[idDEBUG_MEMORY_USE] = 0;
#endif
#define idDEBUG_MEMORY_USE_VERBOSE  7
#ifdef DEBUG_MEMORY_USE_VERBOSE
d_flag[idDEBUG_MEMORY_USE_VERBOSE] = 1;
#else
d_flag[idDEBUG_MEMORY_USE_VERBOSE] = 0;
#endif
#define idDEBUG_TIMESTEP  8
#ifdef DEBUG_TIMESTEP
d_flag[idDEBUG_TIMESTEP] = 1;
#else
d_flag[idDEBUG_TIMESTEP] = 0;
#endif
#define idDEBUG_TIMING  9
#ifdef DEBUG_TIMING
d_flag[idDEBUG_TIMING] = 1;
#else
d_flag[idDEBUG_TIMING] = 0;
#endif
#define idDENSITY_CHUNK_SIZE  10
#ifdef DENSITY_CHUNK_SIZE
d_flag[idDENSITY_CHUNK_SIZE] = 1;
#else
d_flag[idDENSITY_CHUNK_SIZE] = 0;
#endif
#define idELECTRON_ION_NONEQUILIBRIUM  11
#ifdef ELECTRON_ION_NONEQUILIBRIUM
d_flag[idELECTRON_ION_NONEQUILIBRIUM] = 1;
#else
d_flag[idELECTRON_ION_NONEQUILIBRIUM] = 0;
#endif
#define idENRICH  12
#ifdef ENRICH
d_flag[idENRICH] = 1;
#else
d_flag[idENRICH] = 0;
#endif
#define idENRICHMENT  13
#ifdef ENRICHMENT
d_flag[idENRICHMENT] = 1;
#else
d_flag[idENRICHMENT] = 0;
#endif
#define idENRICHMENT_SNIa  14
#ifdef ENRICHMENT_SNIa
d_flag[idENRICHMENT_SNIa] = 1;
#else
d_flag[idENRICHMENT_SNIa] = 0;
#endif
#define idENRICH_SNIa  15
#ifdef ENRICH_SNIa
d_flag[idENRICH_SNIa] = 1;
#else
d_flag[idENRICH_SNIa] = 0;
#endif
#define idFUNCTION  16
#ifdef FUNCTION
d_flag[idFUNCTION] = 1;
#else
d_flag[idFUNCTION] = 0;
#endif
#define idGRAVITY  17
#ifdef GRAVITY
d_flag[idGRAVITY] = 1;
#else
d_flag[idGRAVITY] = 0;
#endif
#define idGRAVITY_IN_RIEMANN  18
#ifdef GRAVITY_IN_RIEMANN
d_flag[idGRAVITY_IN_RIEMANN] = 1;
#else
d_flag[idGRAVITY_IN_RIEMANN] = 0;
#endif
#define idHYDRO  19
#ifdef HYDRO
d_flag[idHYDRO] = 1;
#else
d_flag[idHYDRO] = 0;
#endif
#define idHYDRO_CHUNK_SIZE  20
#ifdef HYDRO_CHUNK_SIZE
d_flag[idHYDRO_CHUNK_SIZE] = 1;
#else
d_flag[idHYDRO_CHUNK_SIZE] = 0;
#endif
#define idHYDRO_TRACERS  21
#ifdef HYDRO_TRACERS
d_flag[idHYDRO_TRACERS] = 1;
#else
d_flag[idHYDRO_TRACERS] = 0;
#endif
#define idHYDRO_TRACERS_NGP  22
#ifdef HYDRO_TRACERS_NGP
d_flag[idHYDRO_TRACERS_NGP] = 1;
#else
d_flag[idHYDRO_TRACERS_NGP] = 0;
#endif
#define idLOG_STAR_CREATION  23
#ifdef LOG_STAR_CREATION
d_flag[idLOG_STAR_CREATION] = 1;
#else
d_flag[idLOG_STAR_CREATION] = 0;
#endif
#define idMOMENTUM_DIFFUSION  24
#ifdef MOMENTUM_DIFFUSION
d_flag[idMOMENTUM_DIFFUSION] = 1;
#else
d_flag[idMOMENTUM_DIFFUSION] = 0;
#endif
#define idMPE_LOG  25
#ifdef MPE_LOG
d_flag[idMPE_LOG] = 1;
#else
d_flag[idMPE_LOG] = 0;
#endif
#define idMPI_MAX_MESSAGE_SIZE  26
#ifdef MPI_MAX_MESSAGE_SIZE
d_flag[idMPI_MAX_MESSAGE_SIZE] = 1;
#else
d_flag[idMPI_MAX_MESSAGE_SIZE] = 0;
#endif
#define idNDEBUG  27
#ifdef NDEBUG
d_flag[idNDEBUG] = 1;
#else
d_flag[idNDEBUG] = 0;
#endif
#define idOLDSTYLE_COOLING_EXPLICIT_SOLVER  28
#ifdef OLDSTYLE_COOLING_EXPLICIT_SOLVER
d_flag[idOLDSTYLE_COOLING_EXPLICIT_SOLVER] = 1;
#else
d_flag[idOLDSTYLE_COOLING_EXPLICIT_SOLVER] = 0;
#endif
#define idOLDSTYLE_PARTICLE_FILE_IGNORE_NGRID  29
#ifdef OLDSTYLE_PARTICLE_FILE_IGNORE_NGRID
d_flag[idOLDSTYLE_PARTICLE_FILE_IGNORE_NGRID] = 1;
#else
d_flag[idOLDSTYLE_PARTICLE_FILE_IGNORE_NGRID] = 0;
#endif
#define idOPENMP_DECLARE_CONST  30
#ifdef OPENMP_DECLARE_CONST
d_flag[idOPENMP_DECLARE_CONST] = 1;
#else
d_flag[idOPENMP_DECLARE_CONST] = 0;
#endif
#define idPARTICLES  31
#ifdef PARTICLES
d_flag[idPARTICLES] = 1;
#else
d_flag[idPARTICLES] = 0;
#endif
#define idPREFIX_JOBNAME_TO_OUTPUT_FILES  32
#ifdef PREFIX_JOBNAME_TO_OUTPUT_FILES
d_flag[idPREFIX_JOBNAME_TO_OUTPUT_FILES] = 1;
#else
d_flag[idPREFIX_JOBNAME_TO_OUTPUT_FILES] = 0;
#endif
#define idRADIATIVE_TRANSFER  33
#ifdef RADIATIVE_TRANSFER
d_flag[idRADIATIVE_TRANSFER] = 1;
#else
d_flag[idRADIATIVE_TRANSFER] = 0;
#endif
#define idREFINEMENT  34
#ifdef REFINEMENT
d_flag[idREFINEMENT] = 1;
#else
d_flag[idREFINEMENT] = 0;
#endif
#define idROOT_GRID_FFT_IJK_TYPE  35
#ifdef ROOT_GRID_FFT_IJK_TYPE
d_flag[idROOT_GRID_FFT_IJK_TYPE] = 1;
#else
d_flag[idROOT_GRID_FFT_IJK_TYPE] = 0;
#endif
#define idRT_1ZONE  36
#ifdef RT_1ZONE
d_flag[idRT_1ZONE] = 1;
#else
d_flag[idRT_1ZONE] = 0;
#endif
#define idRT_8SPECIES  37
#ifdef RT_8SPECIES
d_flag[idRT_8SPECIES] = 1;
#else
d_flag[idRT_8SPECIES] = 0;
#endif
#define idRT_ABSORPTION_CALLBACK_FULL  38
#ifdef RT_ABSORPTION_CALLBACK_FULL
d_flag[idRT_ABSORPTION_CALLBACK_FULL] = 1;
#else
d_flag[idRT_ABSORPTION_CALLBACK_FULL] = 0;
#endif
#define idRT_ADD_EXTERNAL_QSO_BACKGROUND  39
#ifdef RT_ADD_EXTERNAL_QSO_BACKGROUND
d_flag[idRT_ADD_EXTERNAL_QSO_BACKGROUND] = 1;
#else
d_flag[idRT_ADD_EXTERNAL_QSO_BACKGROUND] = 0;
#endif
#define idRT_BACKGROUND_HAARDT_MADAU  40
#ifdef RT_BACKGROUND_HAARDT_MADAU
d_flag[idRT_BACKGROUND_HAARDT_MADAU] = 1;
#else
d_flag[idRT_BACKGROUND_HAARDT_MADAU] = 0;
#endif
#define idRT_BACKGROUND_SELFCONSISTENT  41
#ifdef RT_BACKGROUND_SELFCONSISTENT
d_flag[idRT_BACKGROUND_SELFCONSISTENT] = 1;
#else
d_flag[idRT_BACKGROUND_SELFCONSISTENT] = 0;
#endif
#define idRT_CFI  42
#ifdef RT_CFI
d_flag[idRT_CFI] = 1;
#else
d_flag[idRT_CFI] = 0;
#endif
#define idRT_CHEMISTRY  43
#ifdef RT_CHEMISTRY
d_flag[idRT_CHEMISTRY] = 1;
#else
d_flag[idRT_CHEMISTRY] = 0;
#endif
#define idRT_CHEMISTRY_MINIMAL_MODEL  44
#ifdef RT_CHEMISTRY_MINIMAL_MODEL
d_flag[idRT_CHEMISTRY_MINIMAL_MODEL] = 1;
#else
d_flag[idRT_CHEMISTRY_MINIMAL_MODEL] = 0;
#endif
#define idRT_CUSTOM_DUST_TO_GAS  45
#ifdef RT_CUSTOM_DUST_TO_GAS
d_flag[idRT_CUSTOM_DUST_TO_GAS] = 1;
#else
d_flag[idRT_CUSTOM_DUST_TO_GAS] = 0;
#endif
#define idRT_DEBUG  46
#ifdef RT_DEBUG
d_flag[idRT_DEBUG] = 1;
#else
d_flag[idRT_DEBUG] = 0;
#endif
#define idRT_DEBUG_BLOCK_MASKING  47
#ifdef RT_DEBUG_BLOCK_MASKING
d_flag[idRT_DEBUG_BLOCK_MASKING] = 1;
#else
d_flag[idRT_DEBUG_BLOCK_MASKING] = 0;
#endif
#define idRT_DUST_ABSORBS_BACKGROUND  48
#ifdef RT_DUST_ABSORBS_BACKGROUND
d_flag[idRT_DUST_ABSORBS_BACKGROUND] = 1;
#else
d_flag[idRT_DUST_ABSORBS_BACKGROUND] = 0;
#endif
#define idRT_DUST_CS  49
#ifdef RT_DUST_CS
d_flag[idRT_DUST_CS] = 1;
#else
d_flag[idRT_DUST_CS] = 0;
#endif
#define idRT_EXACT_EOS  50
#ifdef RT_EXACT_EOS
d_flag[idRT_EXACT_EOS] = 1;
#else
d_flag[idRT_EXACT_EOS] = 0;
#endif
#define idRT_EXTERNAL_BACKGROUND  51
#ifdef RT_EXTERNAL_BACKGROUND
d_flag[idRT_EXTERNAL_BACKGROUND] = 1;
#else
d_flag[idRT_EXTERNAL_BACKGROUND] = 0;
#endif
#define idRT_FIXED_ISM  52
#ifdef RT_FIXED_ISM
d_flag[idRT_FIXED_ISM] = 1;
#else
d_flag[idRT_FIXED_ISM] = 0;
#endif
#define idRT_H2_RATE  53
#ifdef RT_H2_RATE
d_flag[idRT_H2_RATE] = 1;
#else
d_flag[idRT_H2_RATE] = 0;
#endif
#define idRT_HIGH_DENSITY  54
#ifdef RT_HIGH_DENSITY
d_flag[idRT_HIGH_DENSITY] = 1;
#else
d_flag[idRT_HIGH_DENSITY] = 0;
#endif
#define idRT_INTERPOLLOG  55
#ifdef RT_INTERPOLLOG
d_flag[idRT_INTERPOLLOG] = 1;
#else
d_flag[idRT_INTERPOLLOG] = 0;
#endif
#define idRT_LWBANDS  56
#ifdef RT_LWBANDS
d_flag[idRT_LWBANDS] = 1;
#else
d_flag[idRT_LWBANDS] = 0;
#endif
#define idRT_LYMAN_ALPHA_HEATING  57
#ifdef RT_LYMAN_ALPHA_HEATING
d_flag[idRT_LYMAN_ALPHA_HEATING] = 1;
#else
d_flag[idRT_LYMAN_ALPHA_HEATING] = 0;
#endif
#define idRT_METHOD_OTVET  58
#ifdef RT_METHOD_OTVET
d_flag[idRT_METHOD_OTVET] = 1;
#else
d_flag[idRT_METHOD_OTVET] = 0;
#endif
#define idRT_NARROW_TABLE  59
#ifdef RT_NARROW_TABLE
d_flag[idRT_NARROW_TABLE] = 1;
#else
d_flag[idRT_NARROW_TABLE] = 0;
#endif
#define idRT_NO_TABLE  60
#ifdef RT_NO_TABLE
d_flag[idRT_NO_TABLE] = 1;
#else
d_flag[idRT_NO_TABLE] = 0;
#endif
#define idRT_OLDSTYLE_H2_SHIELDING  61
#ifdef RT_OLDSTYLE_H2_SHIELDING
d_flag[idRT_OLDSTYLE_H2_SHIELDING] = 1;
#else
d_flag[idRT_OLDSTYLE_H2_SHIELDING] = 0;
#endif
#define idRT_OLDSTYLE_SOURCE_FUNCTION  62
#ifdef RT_OLDSTYLE_SOURCE_FUNCTION
d_flag[idRT_OLDSTYLE_SOURCE_FUNCTION] = 1;
#else
d_flag[idRT_OLDSTYLE_SOURCE_FUNCTION] = 0;
#endif
#define idRT_OTVET_CACHE_ET  63
#ifdef RT_OTVET_CACHE_ET
d_flag[idRT_OTVET_CACHE_ET] = 1;
#else
d_flag[idRT_OTVET_CACHE_ET] = 0;
#endif
#define idRT_OTVET_CACHE_RF  64
#ifdef RT_OTVET_CACHE_RF
d_flag[idRT_OTVET_CACHE_RF] = 1;
#else
d_flag[idRT_OTVET_CACHE_RF] = 0;
#endif
#define idRT_OTVET_SAVE_FLUX  65
#ifdef RT_OTVET_SAVE_FLUX
d_flag[idRT_OTVET_SAVE_FLUX] = 1;
#else
d_flag[idRT_OTVET_SAVE_FLUX] = 0;
#endif
#define idRT_OUTPUT  66
#ifdef RT_OUTPUT
d_flag[idRT_OUTPUT] = 1;
#else
d_flag[idRT_OUTPUT] = 0;
#endif
#define idRT_PAH_CR  67
#ifdef RT_PAH_CR
d_flag[idRT_PAH_CR] = 1;
#else
d_flag[idRT_PAH_CR] = 0;
#endif
#define idRT_PARALLEL_USE_MPI  68
#ifdef RT_PARALLEL_USE_MPI
d_flag[idRT_PARALLEL_USE_MPI] = 1;
#else
d_flag[idRT_PARALLEL_USE_MPI] = 0;
#endif
#define idRT_SINGLE_SOURCE  69
#ifdef RT_SINGLE_SOURCE
d_flag[idRT_SINGLE_SOURCE] = 1;
#else
d_flag[idRT_SINGLE_SOURCE] = 0;
#endif
#define idRT_TRANSFER  70
#ifdef RT_TRANSFER
d_flag[idRT_TRANSFER] = 1;
#else
d_flag[idRT_TRANSFER] = 0;
#endif
#define idRT_TRANSFER_FLUX_CONSERVING  71
#ifdef RT_TRANSFER_FLUX_CONSERVING
d_flag[idRT_TRANSFER_FLUX_CONSERVING] = 1;
#else
d_flag[idRT_TRANSFER_FLUX_CONSERVING] = 0;
#endif
#define idRT_TRANSFER_METHOD  72
#ifdef RT_TRANSFER_METHOD
d_flag[idRT_TRANSFER_METHOD] = 1;
#else
d_flag[idRT_TRANSFER_METHOD] = 0;
#endif
#define idRT_UV  73
#ifdef RT_UV
d_flag[idRT_UV] = 1;
#else
d_flag[idRT_UV] = 0;
#endif
#define idRT_VARIABLE_RF  74
#ifdef RT_VARIABLE_RF
d_flag[idRT_VARIABLE_RF] = 1;
#else
d_flag[idRT_VARIABLE_RF] = 0;
#endif
#define idRT_XLF_BUG_FIX1  75
#ifdef RT_XLF_BUG_FIX1
d_flag[idRT_XLF_BUG_FIX1] = 1;
#else
d_flag[idRT_XLF_BUG_FIX1] = 0;
#endif
#define idRT_XRAYS  76
#ifdef RT_XRAYS
d_flag[idRT_XRAYS] = 1;
#else
d_flag[idRT_XRAYS] = 0;
#endif
#define idSAVE_LOAD_BALANCE_PARTITION  77
#ifdef SAVE_LOAD_BALANCE_PARTITION
d_flag[idSAVE_LOAD_BALANCE_PARTITION] = 1;
#else
d_flag[idSAVE_LOAD_BALANCE_PARTITION] = 0;
#endif
#define idSTARFORM  78
#ifdef STARFORM
d_flag[idSTARFORM] = 1;
#else
d_flag[idSTARFORM] = 0;
#endif
#define idSTAR_FORMATION  79
#ifdef STAR_FORMATION
d_flag[idSTAR_FORMATION] = 1;
#else
d_flag[idSTAR_FORMATION] = 0;
#endif
#define idSTAR_PARTICLE_TYPES  80
#ifdef STAR_PARTICLE_TYPES
d_flag[idSTAR_PARTICLE_TYPES] = 1;
#else
d_flag[idSTAR_PARTICLE_TYPES] = 0;
#endif
#define idUNIQUE_RAND  81
#ifdef UNIQUE_RAND
d_flag[idUNIQUE_RAND] = 1;
#else
d_flag[idUNIQUE_RAND] = 0;
#endif
#define idnum_particle  82
#ifdef num_particle
d_flag[idnum_particle] = 1;
#else
d_flag[idnum_particle] = 0;
#endif
#define PRINT_ALL \
PRINT(AGN) \
PRINT(BLASTWAVE_FEEDBACK) \
PRINT(CONSTANT_TIMESTEP) \
PRINT(COOLING) \
PRINT(COSMOLOGY) \
PRINT(DEBUG) \
PRINT(DEBUG_MEMORY_USE) \
PRINT(DEBUG_MEMORY_USE_VERBOSE) \
PRINT(DEBUG_TIMESTEP) \
PRINT(DEBUG_TIMING) \
PRINT(DENSITY_CHUNK_SIZE) \
PRINT(ELECTRON_ION_NONEQUILIBRIUM) \
PRINT(ENRICH) \
PRINT(ENRICHMENT) \
PRINT(ENRICHMENT_SNIa) \
PRINT(ENRICH_SNIa) \
PRINT(FUNCTION) \
PRINT(GRAVITY) \
PRINT(GRAVITY_IN_RIEMANN) \
PRINT(HYDRO) \
PRINT(HYDRO_CHUNK_SIZE) \
PRINT(HYDRO_TRACERS) \
PRINT(HYDRO_TRACERS_NGP) \
PRINT(LOG_STAR_CREATION) \
PRINT(MOMENTUM_DIFFUSION) \
PRINT(MPE_LOG) \
PRINT(MPI_MAX_MESSAGE_SIZE) \
PRINT(NDEBUG) \
PRINT(OLDSTYLE_COOLING_EXPLICIT_SOLVER) \
PRINT(OLDSTYLE_PARTICLE_FILE_IGNORE_NGRID) \
PRINT(OPENMP_DECLARE_CONST) \
PRINT(PARTICLES) \
PRINT(PREFIX_JOBNAME_TO_OUTPUT_FILES) \
PRINT(RADIATIVE_TRANSFER) \
PRINT(REFINEMENT) \
PRINT(ROOT_GRID_FFT_IJK_TYPE) \
PRINT(RT_1ZONE) \
PRINT(RT_8SPECIES) \
PRINT(RT_ABSORPTION_CALLBACK_FULL) \
PRINT(RT_ADD_EXTERNAL_QSO_BACKGROUND) \
PRINT(RT_BACKGROUND_HAARDT_MADAU) \
PRINT(RT_BACKGROUND_SELFCONSISTENT) \
PRINT(RT_CFI) \
PRINT(RT_CHEMISTRY) \
PRINT(RT_CHEMISTRY_MINIMAL_MODEL) \
PRINT(RT_CUSTOM_DUST_TO_GAS) \
PRINT(RT_DEBUG) \
PRINT(RT_DEBUG_BLOCK_MASKING) \
PRINT(RT_DUST_ABSORBS_BACKGROUND) \
PRINT(RT_DUST_CS) \
PRINT(RT_EXACT_EOS) \
PRINT(RT_EXTERNAL_BACKGROUND) \
PRINT(RT_FIXED_ISM) \
PRINT(RT_H2_RATE) \
PRINT(RT_HIGH_DENSITY) \
PRINT(RT_INTERPOLLOG) \
PRINT(RT_LWBANDS) \
PRINT(RT_LYMAN_ALPHA_HEATING) \
PRINT(RT_METHOD_OTVET) \
PRINT(RT_NARROW_TABLE) \
PRINT(RT_NO_TABLE) \
PRINT(RT_OLDSTYLE_H2_SHIELDING) \
PRINT(RT_OLDSTYLE_SOURCE_FUNCTION) \
PRINT(RT_OTVET_CACHE_ET) \
PRINT(RT_OTVET_CACHE_RF) \
PRINT(RT_OTVET_SAVE_FLUX) \
PRINT(RT_OUTPUT) \
PRINT(RT_PAH_CR) \
PRINT(RT_PARALLEL_USE_MPI) \
PRINT(RT_SINGLE_SOURCE) \
PRINT(RT_TRANSFER) \
PRINT(RT_TRANSFER_FLUX_CONSERVING) \
PRINT(RT_TRANSFER_METHOD) \
PRINT(RT_UV) \
PRINT(RT_VARIABLE_RF) \
PRINT(RT_XLF_BUG_FIX1) \
PRINT(RT_XRAYS) \
PRINT(SAVE_LOAD_BALANCE_PARTITION) \
PRINT(STARFORM) \
PRINT(STAR_FORMATION) \
PRINT(STAR_PARTICLE_TYPES) \
PRINT(UNIQUE_RAND) \
PRINT(num_particle) \
;
