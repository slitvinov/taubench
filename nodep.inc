#ifdef SX3
#pragma vdir nodep
#endif

#ifdef SX4
#pragma vdir nodep
#endif

#ifdef HITACHI
/*poption parallel*/
/*voption indep*/
#endif

#ifdef SUN_SPARC
#ifdef VECTOR
#pragma nomemorydepend
#pragma pipeloop(0)
#endif
#endif
