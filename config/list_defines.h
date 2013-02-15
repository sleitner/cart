#define PRINT(name) \
if(d_flag[id##name]) { fprintf(f,"   %s\n",#name); d_flag[id##name] = 0; }
#define iddef  0
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idENRICHMENT_SNIa  1
#ifdef ENRICHMENT_SNIa
d_flag[idENRICHMENT_SNIa] = 1;
#else
d_flag[idENRICHMENT_SNIa] = 0;
#endif
#define iddef  2
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idPARTICLES  3
#ifdef PARTICLES
d_flag[idPARTICLES] = 1;
#else
d_flag[idPARTICLES] = 0;
#endif
#define iddef  4
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idBLASTWAVE_FEEDBACK  5
#ifdef BLASTWAVE_FEEDBACK
d_flag[idBLASTWAVE_FEEDBACK] = 1;
#else
d_flag[idBLASTWAVE_FEEDBACK] = 0;
#endif
#define iddef  6
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idELECTRON_ION_NONEQUILIBRIUM  7
#ifdef ELECTRON_ION_NONEQUILIBRIUM
d_flag[idELECTRON_ION_NONEQUILIBRIUM] = 1;
#else
d_flag[idELECTRON_ION_NONEQUILIBRIUM] = 0;
#endif
#define iddef  8
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idENRICHMENT  9
#ifdef ENRICHMENT
d_flag[idENRICHMENT] = 1;
#else
d_flag[idENRICHMENT] = 0;
#endif
#define iddef  10
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idFIXED_INTERNAL_ENERGY  11
#ifdef FIXED_INTERNAL_ENERGY
d_flag[idFIXED_INTERNAL_ENERGY] = 1;
#else
d_flag[idFIXED_INTERNAL_ENERGY] = 0;
#endif
#define iddef  12
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idHYDRO  13
#ifdef HYDRO
d_flag[idHYDRO] = 1;
#else
d_flag[idHYDRO] = 0;
#endif
#define iddef  14
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idPARTICLES  15
#ifdef PARTICLES
d_flag[idPARTICLES] = 1;
#else
d_flag[idPARTICLES] = 0;
#endif
#define iddef  16
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idRADIATIVE_TRANSFER  17
#ifdef RADIATIVE_TRANSFER
d_flag[idRADIATIVE_TRANSFER] = 1;
#else
d_flag[idRADIATIVE_TRANSFER] = 0;
#endif
#define iddef  18
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idTURBULENT_ENERGY  19
#ifdef TURBULENT_ENERGY
d_flag[idTURBULENT_ENERGY] = 1;
#else
d_flag[idTURBULENT_ENERGY] = 0;
#endif
#define id!(HYDRO)  20
#ifdef !(HYDRO)
d_flag[id!(HYDRO)] = 1;
#else
d_flag[id!(HYDRO)] = 0;
#endif
#define id&&  21
#ifdef &&
d_flag[id&&] = 1;
#else
d_flag[id&&] = 0;
#endif
#define id(PARTICLES)  22
#ifdef (PARTICLES)
d_flag[id(PARTICLES)] = 1;
#else
d_flag[id(PARTICLES)] = 0;
#endif
#define id(!(PARTICLES)  23
#ifdef (!(PARTICLES)
d_flag[id(!(PARTICLES)] = 1;
#else
d_flag[id(!(PARTICLES)] = 0;
#endif
#define id||  24
#ifdef ||
d_flag[id||] = 1;
#else
d_flag[id||] = 0;
#endif
#define id!(COSMOLOGY))  25
#ifdef !(COSMOLOGY))
d_flag[id!(COSMOLOGY))] = 1;
#else
d_flag[id!(COSMOLOGY))] = 0;
#endif
#define id(COOLING)  26
#ifdef (COOLING)
d_flag[id(COOLING)] = 1;
#else
d_flag[id(COOLING)] = 0;
#endif
#define id&&  27
#ifdef &&
d_flag[id&&] = 1;
#else
d_flag[id&&] = 0;
#endif
#define id!(RADIATIVE_TRANSFER)  28
#ifdef !(RADIATIVE_TRANSFER)
d_flag[id!(RADIATIVE_TRANSFER)] = 1;
#else
d_flag[id!(RADIATIVE_TRANSFER)] = 0;
#endif
#define id(COSMOLOGY)  29
#ifdef (COSMOLOGY)
d_flag[id(COSMOLOGY)] = 1;
#else
d_flag[id(COSMOLOGY)] = 0;
#endif
#define id&&  30
#ifdef &&
d_flag[id&&] = 1;
#else
d_flag[id&&] = 0;
#endif
#define id(PARTICLES)  31
#ifdef (PARTICLES)
d_flag[id(PARTICLES)] = 1;
#else
d_flag[id(PARTICLES)] = 0;
#endif
#define id(ENRICHMENT)  32
#ifdef (ENRICHMENT)
d_flag[id(ENRICHMENT)] = 1;
#else
d_flag[id(ENRICHMENT)] = 0;
#endif
#define id&&  33
#ifdef &&
d_flag[id&&] = 1;
#else
d_flag[id&&] = 0;
#endif
#define id(ENRICHMENT_SNIa)  34
#ifdef (ENRICHMENT_SNIa)
d_flag[id(ENRICHMENT_SNIa)] = 1;
#else
d_flag[id(ENRICHMENT_SNIa)] = 0;
#endif
#define id(GRAVITY)  35
#ifdef (GRAVITY)
d_flag[id(GRAVITY)] = 1;
#else
d_flag[id(GRAVITY)] = 0;
#endif
#define id&&  36
#ifdef &&
d_flag[id&&] = 1;
#else
d_flag[id&&] = 0;
#endif
#define id(!(GRAVITY_IN_RIEMANN))  37
#ifdef (!(GRAVITY_IN_RIEMANN))
d_flag[id(!(GRAVITY_IN_RIEMANN))] = 1;
#else
d_flag[id(!(GRAVITY_IN_RIEMANN))] = 0;
#endif
#define id(GRAVITY)  38
#ifdef (GRAVITY)
d_flag[id(GRAVITY)] = 1;
#else
d_flag[id(GRAVITY)] = 0;
#endif
#define id&&  39
#ifdef &&
d_flag[id&&] = 1;
#else
d_flag[id&&] = 0;
#endif
#define id(PARTICLES)  40
#ifdef (PARTICLES)
d_flag[id(PARTICLES)] = 1;
#else
d_flag[id(PARTICLES)] = 0;
#endif
#define id(GRAVITY)  41
#ifdef (GRAVITY)
d_flag[id(GRAVITY)] = 1;
#else
d_flag[id(GRAVITY)] = 0;
#endif
#define id||  42
#ifdef ||
d_flag[id||] = 1;
#else
d_flag[id||] = 0;
#endif
#define id(RADIATIVE_TRANSFER)  43
#ifdef (RADIATIVE_TRANSFER)
d_flag[id(RADIATIVE_TRANSFER)] = 1;
#else
d_flag[id(RADIATIVE_TRANSFER)] = 0;
#endif
#define id(HYDRO)  44
#ifdef (HYDRO)
d_flag[id(HYDRO)] = 1;
#else
d_flag[id(HYDRO)] = 0;
#endif
#define id&&  45
#ifdef &&
d_flag[id&&] = 1;
#else
d_flag[id&&] = 0;
#endif
#define id(HYDRO_TRACERS)  46
#ifdef (HYDRO_TRACERS)
d_flag[id(HYDRO_TRACERS)] = 1;
#else
d_flag[id(HYDRO_TRACERS)] = 0;
#endif
#define id(HYDRO)  47
#ifdef (HYDRO)
d_flag[id(HYDRO)] = 1;
#else
d_flag[id(HYDRO)] = 0;
#endif
#define id&&  48
#ifdef &&
d_flag[id&&] = 1;
#else
d_flag[id&&] = 0;
#endif
#define id(MOMENTUM_DIFFUSION)  49
#ifdef (MOMENTUM_DIFFUSION)
d_flag[id(MOMENTUM_DIFFUSION)] = 1;
#else
d_flag[id(MOMENTUM_DIFFUSION)] = 0;
#endif
#define id(HYDRO)  50
#ifdef (HYDRO)
d_flag[id(HYDRO)] = 1;
#else
d_flag[id(HYDRO)] = 0;
#endif
#define id&&  51
#ifdef &&
d_flag[id&&] = 1;
#else
d_flag[id&&] = 0;
#endif
#define id(PARTICLES)  52
#ifdef (PARTICLES)
d_flag[id(PARTICLES)] = 1;
#else
d_flag[id(PARTICLES)] = 0;
#endif
#define id(HYDRO)  53
#ifdef (HYDRO)
d_flag[id(HYDRO)] = 1;
#else
d_flag[id(HYDRO)] = 0;
#endif
#define id&&  54
#ifdef &&
d_flag[id&&] = 1;
#else
d_flag[id&&] = 0;
#endif
#define id(STAR_FORMATION)  55
#ifdef (STAR_FORMATION)
d_flag[id(STAR_FORMATION)] = 1;
#else
d_flag[id(STAR_FORMATION)] = 0;
#endif
#define id(HYDRO)  56
#ifdef (HYDRO)
d_flag[id(HYDRO)] = 1;
#else
d_flag[id(HYDRO)] = 0;
#endif
#define id||  57
#ifdef ||
d_flag[id||] = 1;
#else
d_flag[id||] = 0;
#endif
#define id(REFINEMENT)  58
#ifdef (REFINEMENT)
d_flag[id(REFINEMENT)] = 1;
#else
d_flag[id(REFINEMENT)] = 0;
#endif
#define id(PARTICLES)  59
#ifdef (PARTICLES)
d_flag[id(PARTICLES)] = 1;
#else
d_flag[id(PARTICLES)] = 0;
#endif
#define id&&  60
#ifdef &&
d_flag[id&&] = 1;
#else
d_flag[id&&] = 0;
#endif
#define id(STAR_FORMATION)  61
#ifdef (STAR_FORMATION)
d_flag[id(STAR_FORMATION)] = 1;
#else
d_flag[id(STAR_FORMATION)] = 0;
#endif
#define id(RADIATIVE_TRANSFER)  62
#ifdef (RADIATIVE_TRANSFER)
d_flag[id(RADIATIVE_TRANSFER)] = 1;
#else
d_flag[id(RADIATIVE_TRANSFER)] = 0;
#endif
#define id&&  63
#ifdef &&
d_flag[id&&] = 1;
#else
d_flag[id&&] = 0;
#endif
#define id(RT_DEBUG)  64
#ifdef (RT_DEBUG)
d_flag[id(RT_DEBUG)] = 1;
#else
d_flag[id(RT_DEBUG)] = 0;
#endif
#define id(RADIATIVE_TRANSFER)  65
#ifdef (RADIATIVE_TRANSFER)
d_flag[id(RADIATIVE_TRANSFER)] = 1;
#else
d_flag[id(RADIATIVE_TRANSFER)] = 0;
#endif
#define id&&  66
#ifdef &&
d_flag[id&&] = 1;
#else
d_flag[id&&] = 0;
#endif
#define id(RT_TRANSFER)  67
#ifdef (RT_TRANSFER)
d_flag[id(RT_TRANSFER)] = 1;
#else
d_flag[id(RT_TRANSFER)] = 0;
#endif
#define id(RADIATIVE_TRANSFER)  68
#ifdef (RADIATIVE_TRANSFER)
d_flag[id(RADIATIVE_TRANSFER)] = 1;
#else
d_flag[id(RADIATIVE_TRANSFER)] = 0;
#endif
#define id&&  69
#ifdef &&
d_flag[id&&] = 1;
#else
d_flag[id&&] = 0;
#endif
#define id(RT_TRANSFER)  70
#ifdef (RT_TRANSFER)
d_flag[id(RT_TRANSFER)] = 1;
#else
d_flag[id(RT_TRANSFER)] = 0;
#endif
#define id&&  71
#ifdef &&
d_flag[id&&] = 1;
#else
d_flag[id&&] = 0;
#endif
#define id(RT_TRANSFER_METHOD  72
#ifdef (RT_TRANSFER_METHOD
d_flag[id(RT_TRANSFER_METHOD] = 1;
#else
d_flag[id(RT_TRANSFER_METHOD] = 0;
#endif
#define id==  73
#ifdef ==
d_flag[id==] = 1;
#else
d_flag[id==] = 0;
#endif
#define idRT_METHOD_OTVET)  74
#ifdef RT_METHOD_OTVET)
d_flag[idRT_METHOD_OTVET)] = 1;
#else
d_flag[idRT_METHOD_OTVET)] = 0;
#endif
#define id(RT_CFI  75
#ifdef (RT_CFI
d_flag[id(RT_CFI] = 1;
#else
d_flag[id(RT_CFI] = 0;
#endif
#define id==  76
#ifdef ==
d_flag[id==] = 1;
#else
d_flag[id==] = 0;
#endif
#define id0)  77
#ifdef 0)
d_flag[id0)] = 1;
#else
d_flag[id0)] = 0;
#endif
#define id(RT_CFI  78
#ifdef (RT_CFI
d_flag[id(RT_CFI] = 1;
#else
d_flag[id(RT_CFI] = 0;
#endif
#define id==  79
#ifdef ==
d_flag[id==] = 1;
#else
d_flag[id==] = 0;
#endif
#define id1)  80
#ifdef 1)
d_flag[id1)] = 1;
#else
d_flag[id1)] = 0;
#endif
#define id(RT_CFI  81
#ifdef (RT_CFI
d_flag[id(RT_CFI] = 1;
#else
d_flag[id(RT_CFI] = 0;
#endif
#define id==  82
#ifdef ==
d_flag[id==] = 1;
#else
d_flag[id==] = 0;
#endif
#define id2)  83
#ifdef 2)
d_flag[id2)] = 1;
#else
d_flag[id2)] = 0;
#endif
#define id(RT_CFI==1  84
#ifdef (RT_CFI==1
d_flag[id(RT_CFI==1] = 1;
#else
d_flag[id(RT_CFI==1] = 0;
#endif
#define id||  85
#ifdef ||
d_flag[id||] = 1;
#else
d_flag[id||] = 0;
#endif
#define idRT_CFI==2)  86
#ifdef RT_CFI==2)
d_flag[idRT_CFI==2)] = 1;
#else
d_flag[idRT_CFI==2)] = 0;
#endif
#define id(RT_CHEMISTRY)  87
#ifdef (RT_CHEMISTRY)
d_flag[id(RT_CHEMISTRY)] = 1;
#else
d_flag[id(RT_CHEMISTRY)] = 0;
#endif
#define id&&  88
#ifdef &&
d_flag[id&&] = 1;
#else
d_flag[id&&] = 0;
#endif
#define id(RT_LWBANDS)  89
#ifdef (RT_LWBANDS)
d_flag[id(RT_LWBANDS)] = 1;
#else
d_flag[id(RT_LWBANDS)] = 0;
#endif
#define id(RT_DUST_CS  90
#ifdef (RT_DUST_CS
d_flag[id(RT_DUST_CS] = 1;
#else
d_flag[id(RT_DUST_CS] = 0;
#endif
#define id==  91
#ifdef ==
d_flag[id==] = 1;
#else
d_flag[id==] = 0;
#endif
#define id1)  92
#ifdef 1)
d_flag[id1)] = 1;
#else
d_flag[id1)] = 0;
#endif
#define id(RT_EXTERNAL_BACKGROUND  93
#ifdef (RT_EXTERNAL_BACKGROUND
d_flag[id(RT_EXTERNAL_BACKGROUND] = 1;
#else
d_flag[id(RT_EXTERNAL_BACKGROUND] = 0;
#endif
#define id!=  94
#ifdef !=
d_flag[id!=] = 1;
#else
d_flag[id!=] = 0;
#endif
#define idRT_BACKGROUND_SELFCONSISTENT)  95
#ifdef RT_BACKGROUND_SELFCONSISTENT)
d_flag[idRT_BACKGROUND_SELFCONSISTENT)] = 1;
#else
d_flag[idRT_BACKGROUND_SELFCONSISTENT)] = 0;
#endif
#define id(RT_EXTERNAL_BACKGROUND  96
#ifdef (RT_EXTERNAL_BACKGROUND
d_flag[id(RT_EXTERNAL_BACKGROUND] = 1;
#else
d_flag[id(RT_EXTERNAL_BACKGROUND] = 0;
#endif
#define id==  97
#ifdef ==
d_flag[id==] = 1;
#else
d_flag[id==] = 0;
#endif
#define idRT_BACKGROUND_SELFCONSISTENT)  98
#ifdef RT_BACKGROUND_SELFCONSISTENT)
d_flag[idRT_BACKGROUND_SELFCONSISTENT)] = 1;
#else
d_flag[idRT_BACKGROUND_SELFCONSISTENT)] = 0;
#endif
#define id(RT_EXTERNAL_BACKGROUND)  99
#ifdef (RT_EXTERNAL_BACKGROUND)
d_flag[id(RT_EXTERNAL_BACKGROUND)] = 1;
#else
d_flag[id(RT_EXTERNAL_BACKGROUND)] = 0;
#endif
#define id&&  100
#ifdef &&
d_flag[id&&] = 1;
#else
d_flag[id&&] = 0;
#endif
#define id(RT_EXTERNAL_BACKGROUND==RT_BACKGROUND_SELFCONSISTENT)  101
#ifdef (RT_EXTERNAL_BACKGROUND==RT_BACKGROUND_SELFCONSISTENT)
d_flag[id(RT_EXTERNAL_BACKGROUND==RT_BACKGROUND_SELFCONSISTENT)] = 1;
#else
d_flag[id(RT_EXTERNAL_BACKGROUND==RT_BACKGROUND_SELFCONSISTENT)] = 0;
#endif
#define id(RT_EXTERNAL_BACKGROUND==RT_BACKGROUND_HAARDT_MADAU)  102
#ifdef (RT_EXTERNAL_BACKGROUND==RT_BACKGROUND_HAARDT_MADAU)
d_flag[id(RT_EXTERNAL_BACKGROUND==RT_BACKGROUND_HAARDT_MADAU)] = 1;
#else
d_flag[id(RT_EXTERNAL_BACKGROUND==RT_BACKGROUND_HAARDT_MADAU)] = 0;
#endif
#define id(RT_H2_RATE  103
#ifdef (RT_H2_RATE
d_flag[id(RT_H2_RATE] = 1;
#else
d_flag[id(RT_H2_RATE] = 0;
#endif
#define id==  104
#ifdef ==
d_flag[id==] = 1;
#else
d_flag[id==] = 0;
#endif
#define id1)  105
#ifdef 1)
d_flag[id1)] = 1;
#else
d_flag[id1)] = 0;
#endif
#define id(RT_NARROW_TABLE)  106
#ifdef (RT_NARROW_TABLE)
d_flag[id(RT_NARROW_TABLE)] = 1;
#else
d_flag[id(RT_NARROW_TABLE)] = 0;
#endif
#define id||  107
#ifdef ||
d_flag[id||] = 1;
#else
d_flag[id||] = 0;
#endif
#define id(RT_NO_TABLE)  108
#ifdef (RT_NO_TABLE)
d_flag[id(RT_NO_TABLE)] = 1;
#else
d_flag[id(RT_NO_TABLE)] = 0;
#endif
#define id(RT_OUTPUT)  109
#ifdef (RT_OUTPUT)
d_flag[id(RT_OUTPUT)] = 1;
#else
d_flag[id(RT_OUTPUT)] = 0;
#endif
#define id&&  110
#ifdef &&
d_flag[id&&] = 1;
#else
d_flag[id&&] = 0;
#endif
#define id(RT_DEBUG)  111
#ifdef (RT_DEBUG)
d_flag[id(RT_DEBUG)] = 1;
#else
d_flag[id(RT_DEBUG)] = 0;
#endif
#define id(RT_TRANSFER)  112
#ifdef (RT_TRANSFER)
d_flag[id(RT_TRANSFER)] = 1;
#else
d_flag[id(RT_TRANSFER)] = 0;
#endif
#define id(RT_TRANSFER)  113
#ifdef (RT_TRANSFER)
d_flag[id(RT_TRANSFER)] = 1;
#else
d_flag[id(RT_TRANSFER)] = 0;
#endif
#define id&&  114
#ifdef &&
d_flag[id&&] = 1;
#else
d_flag[id&&] = 0;
#endif
#define id(RT_TRANSFER_FLUX_CONSERVING)  115
#ifdef (RT_TRANSFER_FLUX_CONSERVING)
d_flag[id(RT_TRANSFER_FLUX_CONSERVING)] = 1;
#else
d_flag[id(RT_TRANSFER_FLUX_CONSERVING)] = 0;
#endif
#define id(RT_TRANSFER)  116
#ifdef (RT_TRANSFER)
d_flag[id(RT_TRANSFER)] = 1;
#else
d_flag[id(RT_TRANSFER)] = 0;
#endif
#define id&&  117
#ifdef &&
d_flag[id&&] = 1;
#else
d_flag[id&&] = 0;
#endif
#define id(RT_TRANSFER_METHOD  118
#ifdef (RT_TRANSFER_METHOD
d_flag[id(RT_TRANSFER_METHOD] = 1;
#else
d_flag[id(RT_TRANSFER_METHOD] = 0;
#endif
#define id==  119
#ifdef ==
d_flag[id==] = 1;
#else
d_flag[id==] = 0;
#endif
#define idRT_METHOD_OTVET)  120
#ifdef RT_METHOD_OTVET)
d_flag[idRT_METHOD_OTVET)] = 1;
#else
d_flag[idRT_METHOD_OTVET)] = 0;
#endif
#define id(RT_TRANSFER)  121
#ifdef (RT_TRANSFER)
d_flag[id(RT_TRANSFER)] = 1;
#else
d_flag[id(RT_TRANSFER)] = 0;
#endif
#define id&&  122
#ifdef &&
d_flag[id&&] = 1;
#else
d_flag[id&&] = 0;
#endif
#define id(RT_VARIABLE_RF)  123
#ifdef (RT_VARIABLE_RF)
d_flag[id(RT_VARIABLE_RF)] = 1;
#else
d_flag[id(RT_VARIABLE_RF)] = 0;
#endif
#define id(RT_TRANSFER_FLUX_CONSERVING)  124
#ifdef (RT_TRANSFER_FLUX_CONSERVING)
d_flag[id(RT_TRANSFER_FLUX_CONSERVING)] = 1;
#else
d_flag[id(RT_TRANSFER_FLUX_CONSERVING)] = 0;
#endif
#define id||  125
#ifdef ||
d_flag[id||] = 1;
#else
d_flag[id||] = 0;
#endif
#define id(RT_VARIABLE_RF)  126
#ifdef (RT_VARIABLE_RF)
d_flag[id(RT_VARIABLE_RF)] = 1;
#else
d_flag[id(RT_VARIABLE_RF)] = 0;
#endif
#define id(RT_TRANSFER_METHOD  127
#ifdef (RT_TRANSFER_METHOD
d_flag[id(RT_TRANSFER_METHOD] = 1;
#else
d_flag[id(RT_TRANSFER_METHOD] = 0;
#endif
#define id==  128
#ifdef ==
d_flag[id==] = 1;
#else
d_flag[id==] = 0;
#endif
#define idRT_METHOD_OTVET)  129
#ifdef RT_METHOD_OTVET)
d_flag[idRT_METHOD_OTVET)] = 1;
#else
d_flag[idRT_METHOD_OTVET)] = 0;
#endif
#define id(RT_UV)  130
#ifdef (RT_UV)
d_flag[id(RT_UV)] = 1;
#else
d_flag[id(RT_UV)] = 0;
#endif
#define id&&  131
#ifdef &&
d_flag[id&&] = 1;
#else
d_flag[id&&] = 0;
#endif
#define id!(RT_UV_OLDSTYLE_3x1)  132
#ifdef !(RT_UV_OLDSTYLE_3x1)
d_flag[id!(RT_UV_OLDSTYLE_3x1)] = 1;
#else
d_flag[id!(RT_UV_OLDSTYLE_3x1)] = 0;
#endif
#define id(STAR_FORMATION)  133
#ifdef (STAR_FORMATION)
d_flag[id(STAR_FORMATION)] = 1;
#else
d_flag[id(STAR_FORMATION)] = 0;
#endif
#define id&&  134
#ifdef &&
d_flag[id&&] = 1;
#else
d_flag[id&&] = 0;
#endif
#define id(AGN)  135
#ifdef (AGN)
d_flag[id(AGN)] = 1;
#else
d_flag[id(AGN)] = 0;
#endif
#define id(STAR_FORMATION)  136
#ifdef (STAR_FORMATION)
d_flag[id(STAR_FORMATION)] = 1;
#else
d_flag[id(STAR_FORMATION)] = 0;
#endif
#define id&&  137
#ifdef &&
d_flag[id&&] = 1;
#else
d_flag[id&&] = 0;
#endif
#define id(STAR_PARTICLE_TYPES)  138
#ifdef (STAR_PARTICLE_TYPES)
d_flag[id(STAR_PARTICLE_TYPES)] = 1;
#else
d_flag[id(STAR_PARTICLE_TYPES)] = 0;
#endif
#define id(STAR_FORMATION)  139
#ifdef (STAR_FORMATION)
d_flag[id(STAR_FORMATION)] = 1;
#else
d_flag[id(STAR_FORMATION)] = 0;
#endif
#define id&&  140
#ifdef &&
d_flag[id&&] = 1;
#else
d_flag[id&&] = 0;
#endif
#define id(STAR_PARTICLE_TYPES)  141
#ifdef (STAR_PARTICLE_TYPES)
d_flag[id(STAR_PARTICLE_TYPES)] = 1;
#else
d_flag[id(STAR_PARTICLE_TYPES)] = 0;
#endif
#define id(STAR_PARTICLE_TYPES)  142
#ifdef (STAR_PARTICLE_TYPES)
d_flag[id(STAR_PARTICLE_TYPES)] = 1;
#else
d_flag[id(STAR_PARTICLE_TYPES)] = 0;
#endif
#define id&&  143
#ifdef &&
d_flag[id&&] = 1;
#else
d_flag[id&&] = 0;
#endif
#define id(AGN)  144
#ifdef (AGN)
d_flag[id(AGN)] = 1;
#else
d_flag[id(AGN)] = 0;
#endif
#define id(TURBULENT_ENERGY)  145
#ifdef (TURBULENT_ENERGY)
d_flag[id(TURBULENT_ENERGY)] = 1;
#else
d_flag[id(TURBULENT_ENERGY)] = 0;
#endif
#define id||  146
#ifdef ||
d_flag[id||] = 1;
#else
d_flag[id||] = 0;
#endif
#define id(FIXED_INTERNAL_ENERGY)  147
#ifdef (FIXED_INTERNAL_ENERGY)
d_flag[id(FIXED_INTERNAL_ENERGY)] = 1;
#else
d_flag[id(FIXED_INTERNAL_ENERGY)] = 0;
#endif
#define idnum_extra_energy_variables  148
#ifdef num_extra_energy_variables
d_flag[idnum_extra_energy_variables] = 1;
#else
d_flag[idnum_extra_energy_variables] = 0;
#endif
#define id>  149
#ifdef >
d_flag[id>] = 1;
#else
d_flag[id>] = 0;
#endif
#define id0  150
#ifdef 0
d_flag[id0] = 1;
#else
d_flag[id0] = 0;
#endif
#define idnum_extra_energy_variables  151
#ifdef num_extra_energy_variables
d_flag[idnum_extra_energy_variables] = 1;
#else
d_flag[idnum_extra_energy_variables] = 0;
#endif
#define id>  152
#ifdef >
d_flag[id>] = 1;
#else
d_flag[id>] = 0;
#endif
#define id0  153
#ifdef 0
d_flag[id0] = 1;
#else
d_flag[id0] = 0;
#endif
#define iddef  154
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idELECTRON_ION_NONEQUILIBRIUM  155
#ifdef ELECTRON_ION_NONEQUILIBRIUM
d_flag[idELECTRON_ION_NONEQUILIBRIUM] = 1;
#else
d_flag[idELECTRON_ION_NONEQUILIBRIUM] = 0;
#endif
#define iddef  156
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idAGN  157
#ifdef AGN
d_flag[idAGN] = 1;
#else
d_flag[idAGN] = 0;
#endif
#define iddef  158
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idAGN  159
#ifdef AGN
d_flag[idAGN] = 1;
#else
d_flag[idAGN] = 0;
#endif
#define iddef  160
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idAGN  161
#ifdef AGN
d_flag[idAGN] = 1;
#else
d_flag[idAGN] = 0;
#endif
#define iddef  162
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idBLASTWAVE_FEEDBACK  163
#ifdef BLASTWAVE_FEEDBACK
d_flag[idBLASTWAVE_FEEDBACK] = 1;
#else
d_flag[idBLASTWAVE_FEEDBACK] = 0;
#endif
#define iddef  164
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idBLASTWAVE_FEEDBACK  165
#ifdef BLASTWAVE_FEEDBACK
d_flag[idBLASTWAVE_FEEDBACK] = 1;
#else
d_flag[idBLASTWAVE_FEEDBACK] = 0;
#endif
#define iddef  166
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idCONSTANT_TIMESTEP  167
#ifdef CONSTANT_TIMESTEP
d_flag[idCONSTANT_TIMESTEP] = 1;
#else
d_flag[idCONSTANT_TIMESTEP] = 0;
#endif
#define iddef  168
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idCOOLING  169
#ifdef COOLING
d_flag[idCOOLING] = 1;
#else
d_flag[idCOOLING] = 0;
#endif
#define iddef  170
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idCOSMOLOGY  171
#ifdef COSMOLOGY
d_flag[idCOSMOLOGY] = 1;
#else
d_flag[idCOSMOLOGY] = 0;
#endif
#define iddef  172
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idCOSMOLOGY  173
#ifdef COSMOLOGY
d_flag[idCOSMOLOGY] = 1;
#else
d_flag[idCOSMOLOGY] = 0;
#endif
#define iddef  174
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define iddef  175
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idDEBUG_MEMORY_USE  176
#ifdef DEBUG_MEMORY_USE
d_flag[idDEBUG_MEMORY_USE] = 1;
#else
d_flag[idDEBUG_MEMORY_USE] = 0;
#endif
#define iddef  177
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idDEBUG_MEMORY_USE_VERBOSE  178
#ifdef DEBUG_MEMORY_USE_VERBOSE
d_flag[idDEBUG_MEMORY_USE_VERBOSE] = 1;
#else
d_flag[idDEBUG_MEMORY_USE_VERBOSE] = 0;
#endif
#define iddef  179
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define iddef  180
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define iddef  181
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define iddef  182
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define iddef  183
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idELECTRON_ION_NONEQUILIBRIUM  184
#ifdef ELECTRON_ION_NONEQUILIBRIUM
d_flag[idELECTRON_ION_NONEQUILIBRIUM] = 1;
#else
d_flag[idELECTRON_ION_NONEQUILIBRIUM] = 0;
#endif
#define iddef  185
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idELECTRON_ION_NONEQUILIBRIUM  186
#ifdef ELECTRON_ION_NONEQUILIBRIUM
d_flag[idELECTRON_ION_NONEQUILIBRIUM] = 1;
#else
d_flag[idELECTRON_ION_NONEQUILIBRIUM] = 0;
#endif
#define iddef  187
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define iddef  188
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idENRICHMENT  189
#ifdef ENRICHMENT
d_flag[idENRICHMENT] = 1;
#else
d_flag[idENRICHMENT] = 0;
#endif
#define iddef  190
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idENRICHMENT_SNIa  191
#ifdef ENRICHMENT_SNIa
d_flag[idENRICHMENT_SNIa] = 1;
#else
d_flag[idENRICHMENT_SNIa] = 0;
#endif
#define iddef  192
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define iddef  193
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idFIXED_INTERNAL_ENERGY  194
#ifdef FIXED_INTERNAL_ENERGY
d_flag[idFIXED_INTERNAL_ENERGY] = 1;
#else
d_flag[idFIXED_INTERNAL_ENERGY] = 0;
#endif
#define iddef  195
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idGRAVITY  196
#ifdef GRAVITY
d_flag[idGRAVITY] = 1;
#else
d_flag[idGRAVITY] = 0;
#endif
#define iddef  197
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idGRAVITY  198
#ifdef GRAVITY
d_flag[idGRAVITY] = 1;
#else
d_flag[idGRAVITY] = 0;
#endif
#define iddef  199
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idGRAVITY_IN_RIEMANN  200
#ifdef GRAVITY_IN_RIEMANN
d_flag[idGRAVITY_IN_RIEMANN] = 1;
#else
d_flag[idGRAVITY_IN_RIEMANN] = 0;
#endif
#define iddef  201
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idHYDRO  202
#ifdef HYDRO
d_flag[idHYDRO] = 1;
#else
d_flag[idHYDRO] = 0;
#endif
#define iddef  203
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idHYDRO  204
#ifdef HYDRO
d_flag[idHYDRO] = 1;
#else
d_flag[idHYDRO] = 0;
#endif
#define iddef  205
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idHYDRO  206
#ifdef HYDRO
d_flag[idHYDRO] = 1;
#else
d_flag[idHYDRO] = 0;
#endif
#define iddef  207
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idHYDRO  208
#ifdef HYDRO
d_flag[idHYDRO] = 1;
#else
d_flag[idHYDRO] = 0;
#endif
#define iddef  209
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idHYDRO_TRACERS  210
#ifdef HYDRO_TRACERS
d_flag[idHYDRO_TRACERS] = 1;
#else
d_flag[idHYDRO_TRACERS] = 0;
#endif
#define iddef  211
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idHYDRO_TRACERS  212
#ifdef HYDRO_TRACERS
d_flag[idHYDRO_TRACERS] = 1;
#else
d_flag[idHYDRO_TRACERS] = 0;
#endif
#define iddef  213
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idHYDRO_TRACERS_NGP  214
#ifdef HYDRO_TRACERS_NGP
d_flag[idHYDRO_TRACERS_NGP] = 1;
#else
d_flag[idHYDRO_TRACERS_NGP] = 0;
#endif
#define iddef  215
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idLOG_STAR_CREATION  216
#ifdef LOG_STAR_CREATION
d_flag[idLOG_STAR_CREATION] = 1;
#else
d_flag[idLOG_STAR_CREATION] = 0;
#endif
#define iddef  217
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idLOG_STAR_CREATION  218
#ifdef LOG_STAR_CREATION
d_flag[idLOG_STAR_CREATION] = 1;
#else
d_flag[idLOG_STAR_CREATION] = 0;
#endif
#define iddef  219
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idLOG_STAR_CREATION  220
#ifdef LOG_STAR_CREATION
d_flag[idLOG_STAR_CREATION] = 1;
#else
d_flag[idLOG_STAR_CREATION] = 0;
#endif
#define iddef  221
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idMOMENTUM_DIFFUSION  222
#ifdef MOMENTUM_DIFFUSION
d_flag[idMOMENTUM_DIFFUSION] = 1;
#else
d_flag[idMOMENTUM_DIFFUSION] = 0;
#endif
#define iddef  223
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idMOMENTUM_DIFFUSION  224
#ifdef MOMENTUM_DIFFUSION
d_flag[idMOMENTUM_DIFFUSION] = 1;
#else
d_flag[idMOMENTUM_DIFFUSION] = 0;
#endif
#define iddef  225
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idMPE_LOG  226
#ifdef MPE_LOG
d_flag[idMPE_LOG] = 1;
#else
d_flag[idMPE_LOG] = 0;
#endif
#define iddef  227
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idMPI_MAX_MESSAGE_SIZE  228
#ifdef MPI_MAX_MESSAGE_SIZE
d_flag[idMPI_MAX_MESSAGE_SIZE] = 1;
#else
d_flag[idMPI_MAX_MESSAGE_SIZE] = 0;
#endif
#define iddef  229
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idOLDSTYLE_COOLING_EXPLICIT_SOLVER  230
#ifdef OLDSTYLE_COOLING_EXPLICIT_SOLVER
d_flag[idOLDSTYLE_COOLING_EXPLICIT_SOLVER] = 1;
#else
d_flag[idOLDSTYLE_COOLING_EXPLICIT_SOLVER] = 0;
#endif
#define iddef  231
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idOLDSTYLE_SF_ALGORITHM  232
#ifdef OLDSTYLE_SF_ALGORITHM
d_flag[idOLDSTYLE_SF_ALGORITHM] = 1;
#else
d_flag[idOLDSTYLE_SF_ALGORITHM] = 0;
#endif
#define iddef  233
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define iddef  234
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define iddef  235
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idPAPI_PROFILING  236
#ifdef PAPI_PROFILING
d_flag[idPAPI_PROFILING] = 1;
#else
d_flag[idPAPI_PROFILING] = 0;
#endif
#define iddef  237
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idPAPI_PROFILING  238
#ifdef PAPI_PROFILING
d_flag[idPAPI_PROFILING] = 1;
#else
d_flag[idPAPI_PROFILING] = 0;
#endif
#define iddef  239
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idPAPI_PROFILING  240
#ifdef PAPI_PROFILING
d_flag[idPAPI_PROFILING] = 1;
#else
d_flag[idPAPI_PROFILING] = 0;
#endif
#define iddef  241
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idPARTICLES  242
#ifdef PARTICLES
d_flag[idPARTICLES] = 1;
#else
d_flag[idPARTICLES] = 0;
#endif
#define iddef  243
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idPARTICLES  244
#ifdef PARTICLES
d_flag[idPARTICLES] = 1;
#else
d_flag[idPARTICLES] = 0;
#endif
#define iddef  245
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idPARTICLE_VDISP  246
#ifdef PARTICLE_VDISP
d_flag[idPARTICLE_VDISP] = 1;
#else
d_flag[idPARTICLE_VDISP] = 0;
#endif
#define iddef  247
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idPREFIX_JOBNAME_TO_OUTPUT_FILES  248
#ifdef PREFIX_JOBNAME_TO_OUTPUT_FILES
d_flag[idPREFIX_JOBNAME_TO_OUTPUT_FILES] = 1;
#else
d_flag[idPREFIX_JOBNAME_TO_OUTPUT_FILES] = 0;
#endif
#define iddef  249
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idRADIATIVE_TRANSFER  250
#ifdef RADIATIVE_TRANSFER
d_flag[idRADIATIVE_TRANSFER] = 1;
#else
d_flag[idRADIATIVE_TRANSFER] = 0;
#endif
#define iddef  251
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idRADIATIVE_TRANSFER  252
#ifdef RADIATIVE_TRANSFER
d_flag[idRADIATIVE_TRANSFER] = 1;
#else
d_flag[idRADIATIVE_TRANSFER] = 0;
#endif
#define iddef  253
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idRADIATIVE_TRANSFER  254
#ifdef RADIATIVE_TRANSFER
d_flag[idRADIATIVE_TRANSFER] = 1;
#else
d_flag[idRADIATIVE_TRANSFER] = 0;
#endif
#define iddef  255
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idREFINEMENT  256
#ifdef REFINEMENT
d_flag[idREFINEMENT] = 1;
#else
d_flag[idREFINEMENT] = 0;
#endif
#define iddef  257
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define iddef  258
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define iddef  259
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define iddef  260
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idRT_ADD_EXTERNAL_QSO_BACKGROUND  261
#ifdef RT_ADD_EXTERNAL_QSO_BACKGROUND
d_flag[idRT_ADD_EXTERNAL_QSO_BACKGROUND] = 1;
#else
d_flag[idRT_ADD_EXTERNAL_QSO_BACKGROUND] = 0;
#endif
#define iddef  262
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idRT_CHEMISTRY  263
#ifdef RT_CHEMISTRY
d_flag[idRT_CHEMISTRY] = 1;
#else
d_flag[idRT_CHEMISTRY] = 0;
#endif
#define iddef  264
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idRT_CHEMISTRY_MINIMAL_MODEL  265
#ifdef RT_CHEMISTRY_MINIMAL_MODEL
d_flag[idRT_CHEMISTRY_MINIMAL_MODEL] = 1;
#else
d_flag[idRT_CHEMISTRY_MINIMAL_MODEL] = 0;
#endif
#define iddef  266
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idRT_CUSTOM_DUST_TO_GAS  267
#ifdef RT_CUSTOM_DUST_TO_GAS
d_flag[idRT_CUSTOM_DUST_TO_GAS] = 1;
#else
d_flag[idRT_CUSTOM_DUST_TO_GAS] = 0;
#endif
#define iddef  268
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define iddef  269
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idRT_DUST_ABSORBS_BACKGROUND  270
#ifdef RT_DUST_ABSORBS_BACKGROUND
d_flag[idRT_DUST_ABSORBS_BACKGROUND] = 1;
#else
d_flag[idRT_DUST_ABSORBS_BACKGROUND] = 0;
#endif
#define iddef  271
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idRT_EXACT_EOS  272
#ifdef RT_EXACT_EOS
d_flag[idRT_EXACT_EOS] = 1;
#else
d_flag[idRT_EXACT_EOS] = 0;
#endif
#define iddef  273
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idRT_EXTERNAL_BACKGROUND  274
#ifdef RT_EXTERNAL_BACKGROUND
d_flag[idRT_EXTERNAL_BACKGROUND] = 1;
#else
d_flag[idRT_EXTERNAL_BACKGROUND] = 0;
#endif
#define iddef  275
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idRT_FIXED_ISM  276
#ifdef RT_FIXED_ISM
d_flag[idRT_FIXED_ISM] = 1;
#else
d_flag[idRT_FIXED_ISM] = 0;
#endif
#define iddef  277
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idRT_HIGH_DENSITY  278
#ifdef RT_HIGH_DENSITY
d_flag[idRT_HIGH_DENSITY] = 1;
#else
d_flag[idRT_HIGH_DENSITY] = 0;
#endif
#define iddef  279
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define iddef  280
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idRT_LWBANDS  281
#ifdef RT_LWBANDS
d_flag[idRT_LWBANDS] = 1;
#else
d_flag[idRT_LWBANDS] = 0;
#endif
#define iddef  282
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idRT_LYMAN_ALPHA_HEATING  283
#ifdef RT_LYMAN_ALPHA_HEATING
d_flag[idRT_LYMAN_ALPHA_HEATING] = 1;
#else
d_flag[idRT_LYMAN_ALPHA_HEATING] = 0;
#endif
#define iddef  284
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idRT_OLDSTYLE_H2_SHIELDING  285
#ifdef RT_OLDSTYLE_H2_SHIELDING
d_flag[idRT_OLDSTYLE_H2_SHIELDING] = 1;
#else
d_flag[idRT_OLDSTYLE_H2_SHIELDING] = 0;
#endif
#define iddef  286
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idRT_OTVET_CACHE_ET  287
#ifdef RT_OTVET_CACHE_ET
d_flag[idRT_OTVET_CACHE_ET] = 1;
#else
d_flag[idRT_OTVET_CACHE_ET] = 0;
#endif
#define iddef  288
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idRT_OTVET_CACHE_RF  289
#ifdef RT_OTVET_CACHE_RF
d_flag[idRT_OTVET_CACHE_RF] = 1;
#else
d_flag[idRT_OTVET_CACHE_RF] = 0;
#endif
#define iddef  290
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idRT_OTVET_NO_GLOBAL_ARRAY  291
#ifdef RT_OTVET_NO_GLOBAL_ARRAY
d_flag[idRT_OTVET_NO_GLOBAL_ARRAY] = 1;
#else
d_flag[idRT_OTVET_NO_GLOBAL_ARRAY] = 0;
#endif
#define iddef  292
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idRT_OTVET_SAVE_FLUX  293
#ifdef RT_OTVET_SAVE_FLUX
d_flag[idRT_OTVET_SAVE_FLUX] = 1;
#else
d_flag[idRT_OTVET_SAVE_FLUX] = 0;
#endif
#define iddef  294
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idRT_OUTPUT  295
#ifdef RT_OUTPUT
d_flag[idRT_OUTPUT] = 1;
#else
d_flag[idRT_OUTPUT] = 0;
#endif
#define iddef  296
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idRT_PAH_CR  297
#ifdef RT_PAH_CR
d_flag[idRT_PAH_CR] = 1;
#else
d_flag[idRT_PAH_CR] = 0;
#endif
#define iddef  298
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idRT_PARALLEL_USE_MPI  299
#ifdef RT_PARALLEL_USE_MPI
d_flag[idRT_PARALLEL_USE_MPI] = 1;
#else
d_flag[idRT_PARALLEL_USE_MPI] = 0;
#endif
#define iddef  300
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define iddef  301
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idRT_TRANSFER  302
#ifdef RT_TRANSFER
d_flag[idRT_TRANSFER] = 1;
#else
d_flag[idRT_TRANSFER] = 0;
#endif
#define iddef  303
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idRT_TRANSFER  304
#ifdef RT_TRANSFER
d_flag[idRT_TRANSFER] = 1;
#else
d_flag[idRT_TRANSFER] = 0;
#endif
#define iddef  305
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idRT_TRANSFER_FLUX_CONSERVING  306
#ifdef RT_TRANSFER_FLUX_CONSERVING
d_flag[idRT_TRANSFER_FLUX_CONSERVING] = 1;
#else
d_flag[idRT_TRANSFER_FLUX_CONSERVING] = 0;
#endif
#define iddef  307
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idRT_UV  308
#ifdef RT_UV
d_flag[idRT_UV] = 1;
#else
d_flag[idRT_UV] = 0;
#endif
#define iddef  309
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idRT_UV_OLDSTYLE_3x1  310
#ifdef RT_UV_OLDSTYLE_3x1
d_flag[idRT_UV_OLDSTYLE_3x1] = 1;
#else
d_flag[idRT_UV_OLDSTYLE_3x1] = 0;
#endif
#define iddef  311
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define iddef  312
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define iddef  313
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idRT_XRAYS  314
#ifdef RT_XRAYS
d_flag[idRT_XRAYS] = 1;
#else
d_flag[idRT_XRAYS] = 0;
#endif
#define iddef  315
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idSAVE_LOAD_BALANCE_PARTITION  316
#ifdef SAVE_LOAD_BALANCE_PARTITION
d_flag[idSAVE_LOAD_BALANCE_PARTITION] = 1;
#else
d_flag[idSAVE_LOAD_BALANCE_PARTITION] = 0;
#endif
#define iddef  317
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idSF_FEEDBACK  318
#ifdef SF_FEEDBACK
d_flag[idSF_FEEDBACK] = 1;
#else
d_flag[idSF_FEEDBACK] = 0;
#endif
#define iddef  319
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idSF_RECIPE  320
#ifdef SF_RECIPE
d_flag[idSF_RECIPE] = 1;
#else
d_flag[idSF_RECIPE] = 0;
#endif
#define iddef  321
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define iddef  322
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idSTAR_FORMATION  323
#ifdef STAR_FORMATION
d_flag[idSTAR_FORMATION] = 1;
#else
d_flag[idSTAR_FORMATION] = 0;
#endif
#define iddef  324
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idSTAR_FORMATION  325
#ifdef STAR_FORMATION
d_flag[idSTAR_FORMATION] = 1;
#else
d_flag[idSTAR_FORMATION] = 0;
#endif
#define iddef  326
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idSTAR_FORMATION  327
#ifdef STAR_FORMATION
d_flag[idSTAR_FORMATION] = 1;
#else
d_flag[idSTAR_FORMATION] = 0;
#endif
#define iddef  328
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idSTAR_PARTICLE_TYPES  329
#ifdef STAR_PARTICLE_TYPES
d_flag[idSTAR_PARTICLE_TYPES] = 1;
#else
d_flag[idSTAR_PARTICLE_TYPES] = 0;
#endif
#define iddef  330
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idSTAR_PARTICLE_TYPES  331
#ifdef STAR_PARTICLE_TYPES
d_flag[idSTAR_PARTICLE_TYPES] = 1;
#else
d_flag[idSTAR_PARTICLE_TYPES] = 0;
#endif
#define iddef  332
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idTURBULENT_ENERGY  333
#ifdef TURBULENT_ENERGY
d_flag[idTURBULENT_ENERGY] = 1;
#else
d_flag[idTURBULENT_ENERGY] = 0;
#endif
#define iddef  334
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idTURBULENT_ENERGY  335
#ifdef TURBULENT_ENERGY
d_flag[idTURBULENT_ENERGY] = 1;
#else
d_flag[idTURBULENT_ENERGY] = 0;
#endif
#define iddef  336
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idUNIQUE_RAND  337
#ifdef UNIQUE_RAND
d_flag[idUNIQUE_RAND] = 1;
#else
d_flag[idUNIQUE_RAND] = 0;
#endif
#define iddef  338
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idnum_particle  339
#ifdef num_particle
d_flag[idnum_particle] = 1;
#else
d_flag[idnum_particle] = 0;
#endif
#define iddef  340
#ifdef def
d_flag[iddef] = 1;
#else
d_flag[iddef] = 0;
#endif
#define idndef  341
#ifdef ndef
d_flag[idndef] = 1;
#else
d_flag[idndef] = 0;
#endif
#define idCOMPILER_GCC  342
#ifdef COMPILER_GCC
d_flag[idCOMPILER_GCC] = 1;
#else
d_flag[idCOMPILER_GCC] = 0;
#endif
#define idndef  343
#ifdef ndef
d_flag[idndef] = 1;
#else
d_flag[idndef] = 0;
#endif
#define idCOSMOLOGY  344
#ifdef COSMOLOGY
d_flag[idCOSMOLOGY] = 1;
#else
d_flag[idCOSMOLOGY] = 0;
#endif
#define idndef  345
#ifdef ndef
d_flag[idndef] = 1;
#else
d_flag[idndef] = 0;
#endif
#define idDENSITY_CHUNK_SIZE  346
#ifdef DENSITY_CHUNK_SIZE
d_flag[idDENSITY_CHUNK_SIZE] = 1;
#else
d_flag[idDENSITY_CHUNK_SIZE] = 0;
#endif
#define idndef  347
#ifdef ndef
d_flag[idndef] = 1;
#else
d_flag[idndef] = 0;
#endif
#define idENRICHMENT  348
#ifdef ENRICHMENT
d_flag[idENRICHMENT] = 1;
#else
d_flag[idENRICHMENT] = 0;
#endif
#define idndef  349
#ifdef ndef
d_flag[idndef] = 1;
#else
d_flag[idndef] = 0;
#endif
#define idndef  350
#ifdef ndef
d_flag[idndef] = 1;
#else
d_flag[idndef] = 0;
#endif
#define idHYDRO_CHUNK_SIZE  351
#ifdef HYDRO_CHUNK_SIZE
d_flag[idHYDRO_CHUNK_SIZE] = 1;
#else
d_flag[idHYDRO_CHUNK_SIZE] = 0;
#endif
#define idndef  352
#ifdef ndef
d_flag[idndef] = 1;
#else
d_flag[idndef] = 0;
#endif
#define idHYDRO_TRACERS_NGP  353
#ifdef HYDRO_TRACERS_NGP
d_flag[idHYDRO_TRACERS_NGP] = 1;
#else
d_flag[idHYDRO_TRACERS_NGP] = 0;
#endif
#define idndef  354
#ifdef ndef
d_flag[idndef] = 1;
#else
d_flag[idndef] = 0;
#endif
#define idndef  355
#ifdef ndef
d_flag[idndef] = 1;
#else
d_flag[idndef] = 0;
#endif
#define idndef  356
#ifdef ndef
d_flag[idndef] = 1;
#else
d_flag[idndef] = 0;
#endif
#define idM_SQRT_2  357
#ifdef M_SQRT_2
d_flag[idM_SQRT_2] = 1;
#else
d_flag[idM_SQRT_2] = 0;
#endif
#define idndef  358
#ifdef ndef
d_flag[idndef] = 1;
#else
d_flag[idndef] = 0;
#endif
#define idNDEBUG  359
#ifdef NDEBUG
d_flag[idNDEBUG] = 1;
#else
d_flag[idNDEBUG] = 0;
#endif
#define idndef  360
#ifdef ndef
d_flag[idndef] = 1;
#else
d_flag[idndef] = 0;
#endif
#define idOLDSTYLE_COOLING_EXPLICIT_SOLVER  361
#ifdef OLDSTYLE_COOLING_EXPLICIT_SOLVER
d_flag[idOLDSTYLE_COOLING_EXPLICIT_SOLVER] = 1;
#else
d_flag[idOLDSTYLE_COOLING_EXPLICIT_SOLVER] = 0;
#endif
#define idndef  362
#ifdef ndef
d_flag[idndef] = 1;
#else
d_flag[idndef] = 0;
#endif
#define idOLDSTYLE_PARTICLE_FILE_IGNORE_NGRID  363
#ifdef OLDSTYLE_PARTICLE_FILE_IGNORE_NGRID
d_flag[idOLDSTYLE_PARTICLE_FILE_IGNORE_NGRID] = 1;
#else
d_flag[idOLDSTYLE_PARTICLE_FILE_IGNORE_NGRID] = 0;
#endif
#define idndef  364
#ifdef ndef
d_flag[idndef] = 1;
#else
d_flag[idndef] = 0;
#endif
#define idndef  365
#ifdef ndef
d_flag[idndef] = 1;
#else
d_flag[idndef] = 0;
#endif
#define idndef  366
#ifdef ndef
d_flag[idndef] = 1;
#else
d_flag[idndef] = 0;
#endif
#define idndef  367
#ifdef ndef
d_flag[idndef] = 1;
#else
d_flag[idndef] = 0;
#endif
#define idRT_NO_TABLE  368
#ifdef RT_NO_TABLE
d_flag[idRT_NO_TABLE] = 1;
#else
d_flag[idRT_NO_TABLE] = 0;
#endif
#define idndef  369
#ifdef ndef
d_flag[idndef] = 1;
#else
d_flag[idndef] = 0;
#endif
#define PRINT_ALL \
PRINT(def) \
PRINT(ENRICHMENT_SNIa) \
PRINT(def) \
PRINT(PARTICLES) \
PRINT(def) \
PRINT(BLASTWAVE_FEEDBACK) \
PRINT(def) \
PRINT(ELECTRON_ION_NONEQUILIBRIUM) \
PRINT(def) \
PRINT(ENRICHMENT) \
PRINT(def) \
PRINT(FIXED_INTERNAL_ENERGY) \
PRINT(def) \
PRINT(HYDRO) \
PRINT(def) \
PRINT(PARTICLES) \
PRINT(def) \
PRINT(RADIATIVE_TRANSFER) \
PRINT(def) \
PRINT(TURBULENT_ENERGY) \
PRINT(!(HYDRO)) \
PRINT(&&) \
PRINT((PARTICLES)) \
PRINT((!(PARTICLES)) \
PRINT(||) \
PRINT(!(COSMOLOGY))) \
PRINT((COOLING)) \
PRINT(&&) \
PRINT(!(RADIATIVE_TRANSFER)) \
PRINT((COSMOLOGY)) \
PRINT(&&) \
PRINT((PARTICLES)) \
PRINT((ENRICHMENT)) \
PRINT(&&) \
PRINT((ENRICHMENT_SNIa)) \
PRINT((GRAVITY)) \
PRINT(&&) \
PRINT((!(GRAVITY_IN_RIEMANN))) \
PRINT((GRAVITY)) \
PRINT(&&) \
PRINT((PARTICLES)) \
PRINT((GRAVITY)) \
PRINT(||) \
PRINT((RADIATIVE_TRANSFER)) \
PRINT((HYDRO)) \
PRINT(&&) \
PRINT((HYDRO_TRACERS)) \
PRINT((HYDRO)) \
PRINT(&&) \
PRINT((MOMENTUM_DIFFUSION)) \
PRINT((HYDRO)) \
PRINT(&&) \
PRINT((PARTICLES)) \
PRINT((HYDRO)) \
PRINT(&&) \
PRINT((STAR_FORMATION)) \
PRINT((HYDRO)) \
PRINT(||) \
PRINT((REFINEMENT)) \
PRINT((PARTICLES)) \
PRINT(&&) \
PRINT((STAR_FORMATION)) \
PRINT((RADIATIVE_TRANSFER)) \
PRINT(&&) \
PRINT((RT_DEBUG)) \
PRINT((RADIATIVE_TRANSFER)) \
PRINT(&&) \
PRINT((RT_TRANSFER)) \
PRINT((RADIATIVE_TRANSFER)) \
PRINT(&&) \
PRINT((RT_TRANSFER)) \
PRINT(&&) \
PRINT((RT_TRANSFER_METHOD) \
PRINT(==) \
PRINT(RT_METHOD_OTVET)) \
PRINT((RT_CFI) \
PRINT(==) \
PRINT(0)) \
PRINT((RT_CFI) \
PRINT(==) \
PRINT(1)) \
PRINT((RT_CFI) \
PRINT(==) \
PRINT(2)) \
PRINT((RT_CFI==1) \
PRINT(||) \
PRINT(RT_CFI==2)) \
PRINT((RT_CHEMISTRY)) \
PRINT(&&) \
PRINT((RT_LWBANDS)) \
PRINT((RT_DUST_CS) \
PRINT(==) \
PRINT(1)) \
PRINT((RT_EXTERNAL_BACKGROUND) \
PRINT(!=) \
PRINT(RT_BACKGROUND_SELFCONSISTENT)) \
PRINT((RT_EXTERNAL_BACKGROUND) \
PRINT(==) \
PRINT(RT_BACKGROUND_SELFCONSISTENT)) \
PRINT((RT_EXTERNAL_BACKGROUND)) \
PRINT(&&) \
PRINT((RT_EXTERNAL_BACKGROUND==RT_BACKGROUND_SELFCONSISTENT)) \
PRINT((RT_EXTERNAL_BACKGROUND==RT_BACKGROUND_HAARDT_MADAU)) \
PRINT((RT_H2_RATE) \
PRINT(==) \
PRINT(1)) \
PRINT((RT_NARROW_TABLE)) \
PRINT(||) \
PRINT((RT_NO_TABLE)) \
PRINT((RT_OUTPUT)) \
PRINT(&&) \
PRINT((RT_DEBUG)) \
PRINT((RT_TRANSFER)) \
PRINT((RT_TRANSFER)) \
PRINT(&&) \
PRINT((RT_TRANSFER_FLUX_CONSERVING)) \
PRINT((RT_TRANSFER)) \
PRINT(&&) \
PRINT((RT_TRANSFER_METHOD) \
PRINT(==) \
PRINT(RT_METHOD_OTVET)) \
PRINT((RT_TRANSFER)) \
PRINT(&&) \
PRINT((RT_VARIABLE_RF)) \
PRINT((RT_TRANSFER_FLUX_CONSERVING)) \
PRINT(||) \
PRINT((RT_VARIABLE_RF)) \
PRINT((RT_TRANSFER_METHOD) \
PRINT(==) \
PRINT(RT_METHOD_OTVET)) \
PRINT((RT_UV)) \
PRINT(&&) \
PRINT(!(RT_UV_OLDSTYLE_3x1)) \
PRINT((STAR_FORMATION)) \
PRINT(&&) \
PRINT((AGN)) \
PRINT((STAR_FORMATION)) \
PRINT(&&) \
PRINT((STAR_PARTICLE_TYPES)) \
PRINT((STAR_FORMATION)) \
PRINT(&&) \
PRINT((STAR_PARTICLE_TYPES)) \
PRINT((STAR_PARTICLE_TYPES)) \
PRINT(&&) \
PRINT((AGN)) \
PRINT((TURBULENT_ENERGY)) \
PRINT(||) \
PRINT((FIXED_INTERNAL_ENERGY)) \
PRINT(num_extra_energy_variables) \
PRINT(>) \
PRINT(0) \
PRINT(num_extra_energy_variables) \
PRINT(>) \
PRINT(0) \
PRINT(def) \
PRINT(ELECTRON_ION_NONEQUILIBRIUM) \
PRINT(def) \
PRINT(AGN) \
PRINT(def) \
PRINT(AGN) \
PRINT(def) \
PRINT(AGN) \
PRINT(def) \
PRINT(BLASTWAVE_FEEDBACK) \
PRINT(def) \
PRINT(BLASTWAVE_FEEDBACK) \
PRINT(def) \
PRINT(CONSTANT_TIMESTEP) \
PRINT(def) \
PRINT(COOLING) \
PRINT(def) \
PRINT(COSMOLOGY) \
PRINT(def) \
PRINT(COSMOLOGY) \
PRINT(def) \
PRINT(def) \
PRINT(DEBUG_MEMORY_USE) \
PRINT(def) \
PRINT(DEBUG_MEMORY_USE_VERBOSE) \
PRINT(def) \
PRINT(def) \
PRINT(def) \
PRINT(def) \
PRINT(def) \
PRINT(ELECTRON_ION_NONEQUILIBRIUM) \
PRINT(def) \
PRINT(ELECTRON_ION_NONEQUILIBRIUM) \
PRINT(def) \
PRINT(def) \
PRINT(ENRICHMENT) \
PRINT(def) \
PRINT(ENRICHMENT_SNIa) \
PRINT(def) \
PRINT(def) \
PRINT(FIXED_INTERNAL_ENERGY) \
PRINT(def) \
PRINT(GRAVITY) \
PRINT(def) \
PRINT(GRAVITY) \
PRINT(def) \
PRINT(GRAVITY_IN_RIEMANN) \
PRINT(def) \
PRINT(HYDRO) \
PRINT(def) \
PRINT(HYDRO) \
PRINT(def) \
PRINT(HYDRO) \
PRINT(def) \
PRINT(HYDRO) \
PRINT(def) \
PRINT(HYDRO_TRACERS) \
PRINT(def) \
PRINT(HYDRO_TRACERS) \
PRINT(def) \
PRINT(HYDRO_TRACERS_NGP) \
PRINT(def) \
PRINT(LOG_STAR_CREATION) \
PRINT(def) \
PRINT(LOG_STAR_CREATION) \
PRINT(def) \
PRINT(LOG_STAR_CREATION) \
PRINT(def) \
PRINT(MOMENTUM_DIFFUSION) \
PRINT(def) \
PRINT(MOMENTUM_DIFFUSION) \
PRINT(def) \
PRINT(MPE_LOG) \
PRINT(def) \
PRINT(MPI_MAX_MESSAGE_SIZE) \
PRINT(def) \
PRINT(OLDSTYLE_COOLING_EXPLICIT_SOLVER) \
PRINT(def) \
PRINT(OLDSTYLE_SF_ALGORITHM) \
PRINT(def) \
PRINT(def) \
PRINT(def) \
PRINT(PAPI_PROFILING) \
PRINT(def) \
PRINT(PAPI_PROFILING) \
PRINT(def) \
PRINT(PAPI_PROFILING) \
PRINT(def) \
PRINT(PARTICLES) \
PRINT(def) \
PRINT(PARTICLES) \
PRINT(def) \
PRINT(PARTICLE_VDISP) \
PRINT(def) \
PRINT(PREFIX_JOBNAME_TO_OUTPUT_FILES) \
PRINT(def) \
PRINT(RADIATIVE_TRANSFER) \
PRINT(def) \
PRINT(RADIATIVE_TRANSFER) \
PRINT(def) \
PRINT(RADIATIVE_TRANSFER) \
PRINT(def) \
PRINT(REFINEMENT) \
PRINT(def) \
PRINT(def) \
PRINT(def) \
PRINT(def) \
PRINT(RT_ADD_EXTERNAL_QSO_BACKGROUND) \
PRINT(def) \
PRINT(RT_CHEMISTRY) \
PRINT(def) \
PRINT(RT_CHEMISTRY_MINIMAL_MODEL) \
PRINT(def) \
PRINT(RT_CUSTOM_DUST_TO_GAS) \
PRINT(def) \
PRINT(def) \
PRINT(RT_DUST_ABSORBS_BACKGROUND) \
PRINT(def) \
PRINT(RT_EXACT_EOS) \
PRINT(def) \
PRINT(RT_EXTERNAL_BACKGROUND) \
PRINT(def) \
PRINT(RT_FIXED_ISM) \
PRINT(def) \
PRINT(RT_HIGH_DENSITY) \
PRINT(def) \
PRINT(def) \
PRINT(RT_LWBANDS) \
PRINT(def) \
PRINT(RT_LYMAN_ALPHA_HEATING) \
PRINT(def) \
PRINT(RT_OLDSTYLE_H2_SHIELDING) \
PRINT(def) \
PRINT(RT_OTVET_CACHE_ET) \
PRINT(def) \
PRINT(RT_OTVET_CACHE_RF) \
PRINT(def) \
PRINT(RT_OTVET_NO_GLOBAL_ARRAY) \
PRINT(def) \
PRINT(RT_OTVET_SAVE_FLUX) \
PRINT(def) \
PRINT(RT_OUTPUT) \
PRINT(def) \
PRINT(RT_PAH_CR) \
PRINT(def) \
PRINT(RT_PARALLEL_USE_MPI) \
PRINT(def) \
PRINT(def) \
PRINT(RT_TRANSFER) \
PRINT(def) \
PRINT(RT_TRANSFER) \
PRINT(def) \
PRINT(RT_TRANSFER_FLUX_CONSERVING) \
PRINT(def) \
PRINT(RT_UV) \
PRINT(def) \
PRINT(RT_UV_OLDSTYLE_3x1) \
PRINT(def) \
PRINT(def) \
PRINT(def) \
PRINT(RT_XRAYS) \
PRINT(def) \
PRINT(SAVE_LOAD_BALANCE_PARTITION) \
PRINT(def) \
PRINT(SF_FEEDBACK) \
PRINT(def) \
PRINT(SF_RECIPE) \
PRINT(def) \
PRINT(def) \
PRINT(STAR_FORMATION) \
PRINT(def) \
PRINT(STAR_FORMATION) \
PRINT(def) \
PRINT(STAR_FORMATION) \
PRINT(def) \
PRINT(STAR_PARTICLE_TYPES) \
PRINT(def) \
PRINT(STAR_PARTICLE_TYPES) \
PRINT(def) \
PRINT(TURBULENT_ENERGY) \
PRINT(def) \
PRINT(TURBULENT_ENERGY) \
PRINT(def) \
PRINT(UNIQUE_RAND) \
PRINT(def) \
PRINT(num_particle) \
PRINT(def) \
PRINT(ndef) \
PRINT(COMPILER_GCC) \
PRINT(ndef) \
PRINT(COSMOLOGY) \
PRINT(ndef) \
PRINT(DENSITY_CHUNK_SIZE) \
PRINT(ndef) \
PRINT(ENRICHMENT) \
PRINT(ndef) \
PRINT(ndef) \
PRINT(HYDRO_CHUNK_SIZE) \
PRINT(ndef) \
PRINT(HYDRO_TRACERS_NGP) \
PRINT(ndef) \
PRINT(ndef) \
PRINT(ndef) \
PRINT(M_SQRT_2) \
PRINT(ndef) \
PRINT(NDEBUG) \
PRINT(ndef) \
PRINT(OLDSTYLE_COOLING_EXPLICIT_SOLVER) \
PRINT(ndef) \
PRINT(OLDSTYLE_PARTICLE_FILE_IGNORE_NGRID) \
PRINT(ndef) \
PRINT(ndef) \
PRINT(ndef) \
PRINT(ndef) \
PRINT(RT_NO_TABLE) \
PRINT(ndef) \
;
