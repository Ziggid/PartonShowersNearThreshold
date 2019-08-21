/* -------------------------------------------------------------------------  
 * Name            : rng.h  (header file for the library file rng.c) 
 * Author          : Steve Park & Dave Geyer  
 * Language        : ANSI C 
 * Latest Revision : 09-11-98
 * ------------------------------------------------------------------------- 
 */

#ifndef _RNG_
#define _RNG_

#define MODULUS    2147483647L /* DON'T CHANGE THIS VALUE                   */
#define MULTIPLIER 48271L      /* use 16807 for the "minimal standard"      */
#define CHECK      399268537L  /* use 1043616065 for the "minimal standard" */
#define DEFAULT    123456789L  /* initial seed, use 0 < DEFAULT < MODULUS   */

double Random(void);
void   GetSeed(long *x);
void   PutSeed(long x);
void   TestRandom(void);

#endif //_RNG_
