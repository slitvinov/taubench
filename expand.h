#ifndef EXPAND_COUNT
#error 'EXPAND_COUNT' is not defined
#endif

#ifdef SX3
#pragma vdir expand = EXPAND_COUNT
#endif

#ifdef SX4

#if EXPAND_COUNT == 0
#define _UNROLL_ 0
#endif
#if EXPAND_COUNT == 1
#define _UNROLL_ 1
#endif
#if EXPAND_COUNT == 2
#define _UNROLL_ 2
#endif
#if EXPAND_COUNT == 3
#define _UNROLL_ 3
#endif
#if EXPAND_COUNT == 4
#define _UNROLL_ 4
#endif
#if EXPAND_COUNT == 5
#define _UNROLL_ 5
#endif
#if EXPAND_COUNT == 6
#define _UNROLL_ 6
#endif
#if EXPAND_COUNT == 7
#define _UNROLL_ 7
#endif
#if EXPAND_COUNT == 8
#define _UNROLL_ 8
#endif
#if EXPAND_COUNT == 9
#define _UNROLL_ 9
#endif
#if EXPAND_COUNT == 10
#define _UNROLL_ 10
#endif
#if EXPAND_COUNT > 10
#define _UNROLL_ EXPAND_COUNT
#endif

#pragma vdir expand = _UNROLL_
#undef _UNROLL_
#endif

#ifdef HITACHI
#if EXPAND_COUNT == 1
/*soption unroll(1) */
#endif
#if EXPAND_COUNT == 2
/*soption unroll(2) */
#endif
#if EXPAND_COUNT == 3
/*soption unroll(3) */
#endif
#if EXPAND_COUNT == 4
/*soption unroll(4) */
#endif
#if EXPAND_COUNT == 5
/*soption unroll(5) */
#endif
#if EXPAND_COUNT == 6
/*soption unroll(6) */
#endif
#if EXPAND_COUNT == 7
/*soption unroll(7) */
#endif
#if EXPAND_COUNT == 8
/*soption unroll(8) */
#endif
#if EXPAND_COUNT == 9
/*soption unroll(9) */
#endif
#if EXPAND_COUNT == 10
/*soption unroll(10) */
#endif
#if EXPAND_COUNT == 11
/*soption unroll(11) */
#endif
#if EXPAND_COUNT == 12
/*soption unroll(12) */
#endif
#if EXPAND_COUNT == 13
/*soption unroll(13) */
#endif
#if EXPAND_COUNT == 14
/*soption unroll(14) */
#endif
#if EXPAND_COUNT == 15
/*soption unroll(15) */
#endif
#if EXPAND_COUNT == 16
/*soption unroll(16) */
#endif
#endif

#ifdef SUN_SPARC

#if EXPAND_COUNT == 0
#define _UNROLL_ 0
#endif
#if EXPAND_COUNT == 1
#define _UNROLL_ 1
#endif
#if EXPAND_COUNT == 2
#define _UNROLL_ 2
#endif
#if EXPAND_COUNT == 3
#define _UNROLL_ 3
#endif
#if EXPAND_COUNT == 4
#define _UNROLL_ 4
#endif
#if EXPAND_COUNT == 5
#define _UNROLL_ 5
#endif
#if EXPAND_COUNT == 6
#define _UNROLL_ 6
#endif
#if EXPAND_COUNT == 7
#define _UNROLL_ 7
#endif
#if EXPAND_COUNT == 8
#define _UNROLL_ 8
#endif
#if EXPAND_COUNT == 9
#define _UNROLL_ 9
#endif
#if EXPAND_COUNT == 10
#define _UNROLL_ 10
#endif
#if EXPAND_COUNT > 10
#define _UNROLL_ EXPAND_COUNT
#endif

#if _UNROLL_ > 0
#pragma unroll(_UNROLL_)
#endif
#undef _UNROLL_
#endif

#undef EXPAND_COUNT
