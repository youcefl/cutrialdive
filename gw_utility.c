/*
Nomenclature:
m mpz
w woltman
g giants
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <gmp.h>
#include "giants.h"
#include "gwnum.h"
#include "gwthread.h"
#include "gwcommon.h"
#include "gw_utility.h"

#define MAX_STRING_LENGTH 100000
#define GW_THRESHOLD 4096
#define SAFE 1024 // Many careful loops for large base loop.
#define STEP 256 // Iterations before check.
#define m2w(mx,gx,wx,string) mpz_get_str(string,10,mx);ctog(string,gx);gianttogw(gwdata,gx,wx);
#define w2m(wx,gx,mx) gwtogiant(gwdata,wx,gx);gtompz(gx,mx);

/*****************************************************************************/

void small_base_exp ( mpz_ptr mr, double db, mpz_ptr me, mpz_srcptr mn )
{
    int fft_size = 1;
    int g_length;
    int j;
    char string [MAX_STRING_LENGTH];
    gwnum wr, wc;
    gwhandle *gwdata;
    giant gr, gn;
    mpz_get_str ( string, 10, mn );
    g_length = ( strlen ( string ) >> 2 ) + 8;
    gr  = newgiant ( g_length );
    gn  = newgiant ( g_length );
    ctog( string, gn );
    gwdata = (gwhandle*) malloc ( sizeof ( gwhandle ) );
    gwinit ( gwdata );
    gwset_maxmulbyconst ( gwdata, (long) db );
    gwset_larger_fftlen_count ( gwdata, fft_size );
    gwsetup_general_mod_giant ( gwdata, gn );
    gwsetmulbyconst ( gwdata, (long) db );
    wr = gwalloc ( gwdata );
    wc = gwalloc ( gwdata );
    gwerror_checking ( gwdata, 1 );
    mpz_set_ui ( mr, (unsigned) db );
    m2w ( mr, gr, wr, string );
    for ( int i = mpz_sizeinbase ( me, 2 ) - 2; i >= 0; )
    {
        gwcopy ( gwdata, wr, wc );
        for ( j = 1; j <= STEP && i >= 0; j++, i-- )
        {
            if ( mpz_tstbit ( me, i ) )
            {
                gwsquare2 ( gwdata, wr, wr, 0x0200 );
            }
            else
            {
                gwsquare2 ( gwdata, wr, wr, 0x0000 );
            }
        }
        if ( gw_get_maxerr ( gwdata ) > 0.44 )
        {
            printf ( "*** (sm) Excessive round off error %lf. Increasing FFT size at bit %d.\n", gw_get_maxerr( gwdata ), i );
            gwtogiant ( gwdata, wc, gr );
            gwfreeall ( gwdata );
            gwinit ( gwdata );
            gwset_maxmulbyconst ( gwdata, (long) db );
            gwset_larger_fftlen_count ( gwdata, ++fft_size );
            gwsetup_general_mod_giant ( gwdata, gn );
            gwsetmulbyconst ( gwdata, (long) db );
            wr = gwalloc ( gwdata );
            wc = gwalloc ( gwdata );
            gwerror_checking ( gwdata, 1 );
            gianttogw ( gwdata, gr, wr );
            i += j;
        }
    }  
    w2m ( wr, gr, mr );
    gwdone ( gwdata );
}

/*****************************************************************************/

void medium_base_exp ( mpz_ptr mr, double db, mpz_ptr me, mpz_srcptr mn )
{
    int fft_size = 1;
    int g_length;
    int j;
    char string [MAX_STRING_LENGTH];
    gwnum wr, wc;
    gwhandle *gwdata;
    giant gr, gn;
    mpz_get_str ( string, 10, mn );
    g_length = ( strlen ( string ) >> 2 ) + 8;
    gr = newgiant ( g_length );
    gn = newgiant ( g_length );
    ctog( string, gn );
    gwdata = (gwhandle*) malloc ( sizeof ( gwhandle ) );
    gwinit ( gwdata );
    gwset_larger_fftlen_count ( gwdata, fft_size );
    gwsetup_general_mod_giant ( gwdata, gn );
    wr = gwalloc ( gwdata );
    wc = gwalloc ( gwdata );
    gwerror_checking ( gwdata, 1 );
    mpz_set_ui ( mr, (unsigned) db );
    m2w ( mr, gr, wr, string );
    for ( int i = mpz_sizeinbase ( me, 2 ) - 2 ; i >= 0; )
    {
        gwcopy ( gwdata, wr, wc );
        for ( j = 1; j <= STEP && i >= 0; j++, i-- )
        {
            gwsquare2 ( gwdata, wr, wr, 0x0000 );
            if ( mpz_tstbit ( me, i ) )
            {
                gwsmallmul ( gwdata, db, wr );
            }
        }
        if ( gw_get_maxerr ( gwdata ) > 0.44 )
        {
            printf ( "*** (md) Excessive round off error %lf. Increasing FFT size at bit %d.\n", gw_get_maxerr( gwdata ), i );
            gwtogiant ( gwdata, wc, gr );
            gwfreeall ( gwdata );
            gwinit ( gwdata );
            gwset_larger_fftlen_count ( gwdata, ++fft_size );
            gwsetup_general_mod_giant ( gwdata, gn );
            wr = gwalloc ( gwdata );
            wc = gwalloc ( gwdata );
            gwerror_checking ( gwdata, 1 );
            gianttogw ( gwdata, gr, wr );
            i += j;
        }
    }
    w2m ( wr, gr, mr );
    gwdone ( gwdata );
}

/*****************************************************************************/

void large_base_exp ( mpz_ptr mr, mpz_ptr mb, mpz_ptr me, mpz_srcptr mn )
{
    int fft_size = 1;
    int g_length;
    int j = mpz_sizeinbase ( me, 2 ) - 2;
    int k = j - SAFE;
    char string [MAX_STRING_LENGTH];
    giant gr, gb, gn;
    mpz_get_str ( string, 10, mn );
    g_length = ( strlen ( string ) >> 2 ) + 8;
    gr = newgiant ( g_length );
    gb = newgiant ( g_length );
    gn = newgiant ( g_length );
    ctog ( string, gn );
    gwnum wr, wb, wc;
    gwhandle *gwdata;
    gwdata = (gwhandle*) malloc ( sizeof ( gwhandle ) );
    gwinit ( gwdata );
    gwset_larger_fftlen_count ( gwdata, fft_size );
    gwsetup_general_mod_giant ( gwdata, gn );
    wr = gwalloc ( gwdata );
    wb = gwalloc ( gwdata );
    wc = gwalloc ( gwdata );
    gwerror_checking ( gwdata, 1 );
    m2w ( mb, gb, wb, string );
    gwcopy ( gwdata, wb, wr );
    for ( int i = j; i > k; i-- )
    {
        gwmul3_carefully ( gwdata, wr, wr, wr, 0x0000 );
        if ( mpz_tstbit ( me, i ) )
        {
            gwmul_carefully ( gwdata, wb, wr );
        }
    }
    for ( int i = k ; i > SAFE; )
    {
        gwcopy ( gwdata, wr, wc );
        for ( j = 1; j <= STEP && i > SAFE; j++, i-- )
        {
            gwsquare2 ( gwdata, wr, wr, 0x0000 );
            if ( mpz_tstbit ( me, i ) )
            {
                gwmul3 ( gwdata, wb, wr, wr, 0x0000 );
            }
        }
        if ( gw_get_maxerr ( gwdata ) > 0.44 )
        {
            printf ( "*** (lg) Excessive round off error %lf. Increasing FFT size at bit %d.\n", gw_get_maxerr( gwdata ), i );
            gwtogiant ( gwdata, wc, gr );
            gwfreeall ( gwdata );
            gwinit ( gwdata );
            gwset_larger_fftlen_count ( gwdata, ++fft_size );
            gwsetup_general_mod_giant ( gwdata, gn );
            wr = gwalloc ( gwdata );
            wb = gwalloc ( gwdata );
            wc = gwalloc ( gwdata );
            gwerror_checking ( gwdata, 1 );
            gianttogw ( gwdata, gb, wb );
            gianttogw ( gwdata, gr, wr );
            i += j;
        }
    }
    for ( int i = SAFE; i >= 0; i-- )
    {
        gwmul3_carefully ( gwdata, wr, wr, wr, 0x0000 );
        if ( mpz_tstbit ( me, i ) )
        {
            gwmul_carefully ( gwdata, wb, wr );
        }
    }
    w2m ( wr, gr, mr );
    gwdone ( gwdata );
}

/*****************************************************************************/

int gw_prp ( mpz_t mn )
{ // base 3 euler
    if ( mpz_sizeinbase ( mn, 2 ) < GW_THRESHOLD )
    {
        return ( mpz_probab_prime_p ( mn, 0) > 0 );
    }
    int res = 0;
    mpz_t mr;
    mpz_t me;
    mpz_init ( mr );
    mpz_init ( me );
    mpz_sub_ui ( me, mn, 1 );
    mpz_tdiv_q_2exp ( me, me, 1 );
    small_base_exp ( mr, 3.0, me, mn );
    if ( mpz_cmp_ui ( mr, 1 ) == 0 )
    {
        res = 1;
    }
    else
    {
        mpz_add_ui ( mr, mr, 1 );
        mpz_mod (mr, mr, mn );
        if ( mpz_cmp_ui ( mr, 0 ) == 0 )
        {
            res = 1;
        }
    }
    mpz_clear ( mr );
    mpz_clear ( me );
    return ( res );
}

/*****************************************************************************/

void gw_powm ( mpz_ptr mr, mpz_ptr mb_inp, mpz_ptr me, mpz_srcptr mn )
{    
    if ( mpz_sizeinbase ( me, 2 ) < GW_THRESHOLD )
    {
        mpz_powm ( mr, mb_inp, me, mn );
        return;
    }
    int base_sign = mpz_sgn ( mb_inp );
    int base_size = mpz_sizeinbase ( mb_inp, 2 );
    double db;
    mpz_t mb;
    mpz_init_set ( mb, mb_inp );
    if ( base_sign == -1 )
    {
        mpz_neg ( mb, mb );
    }
    if ( base_size < 26 )
    {
        db = (double) mpz_get_ui ( mb );
        if ( base_size < 9 )
        {
            small_base_exp ( mr, db, me, mn );
        }
        else
        {

            medium_base_exp ( mr, db, me, mn );
        }
    }
    else
    {
        large_base_exp ( mr, mb, me, mn );
    }
    if ( base_sign == -1 && mpz_odd_p ( me ) )
    {
        mpz_neg ( mr, mr );
        mpz_mod ( mr, mr, mn );
    }
    mpz_clear ( mb );
}

/*****************************************************************************/
/*****************************************************************************/
