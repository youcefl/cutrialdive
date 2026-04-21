/*
Nomenclature:
m mpz
w woltman
g giants
*/

#include <gmp.h>

#ifdef __cplusplus
extern "C" {
#endif

void small_base_exp ( mpz_ptr mr, double db, mpz_ptr me, mpz_srcptr mn );
void medium_base_exp ( mpz_ptr mr, double db, mpz_ptr me, mpz_srcptr mn );
void large_base_exp ( mpz_ptr mr, mpz_ptr mb, mpz_ptr me, mpz_srcptr mn );
int gw_prp ( mpz_t mn );
void gw_powm ( mpz_ptr mr, mpz_ptr mb_inp, mpz_ptr me, mpz_srcptr mn );

/*****************************************************************************/
/*****************************************************************************/

#ifdef __cplusplus
};
#endif
