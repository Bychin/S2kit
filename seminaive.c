/***************************************************************************
  **************************************************************************
  
                           S2kit 1.0

          A lite version of Spherical Harmonic Transform Kit

   Peter Kostelec, Dan Rockmore
   {geelong,rockmore}@cs.dartmouth.edu
  
   Contact: Peter Kostelec
            geelong@cs.dartmouth.edu
  
   Copyright 2004 Peter Kostelec, Dan Rockmore

   This file is part of S2kit.

   S2kit is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   S2kit is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with S2kit; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

   See the accompanying LICENSE file for details.
  
  ************************************************************************
  ************************************************************************/


/* Source code for computing the Legendre transform where
   projections are carried out in cosine space, i.e., the
   "seminaive" algorithm.

   For a description, see the related paper or Sean's thesis.

*/

#include <math.h>
#include <stdio.h>
#include <string.h>   /** for memcpy **/
#include "fftw3.h"

#include "cospmls.h"


/************************************************************************/
/* InvSemiNaiveReduced computes the inverse Legendre transform
   using the transposed seminaive algorithm.  Note that because
   the Legendre transform is orthogonal, the inverse can be
   computed by transposing the matrix formulation of the
   problem.

   The forward transform looks like

   l = PCWf

   where f is the data vector, W is a quadrature matrix,
   C is a cosine transform matrix, P is a matrix
   full of coefficients of the cosine series representation
   of each Pml function P(m,m) P(m,m+1) ... P(m,bw-1),
   and l is the (associated) Legendre series representation
   of f.

   So to do the inverse, you do

   f = trans(C) trans(P) l

   so you need to transpose the matrix P from the forward transform
   and then do a cosine series evaluation.  No quadrature matrix
   is necessary.  If order m is odd, then there is also a sin
   factor that needs to be accounted for.

   Note that this function was written to be part of a full
   spherical harmonic transform, so a lot of precomputation
   has been assumed.

   Input argument description

   coeffs - a double pointer to an array of length
            (bw-m) containing associated
	    Legendre series coefficients.  Assumed
	    that first entry contains the P(m,m)
	    coefficient.

   bw - problem bandwidth

   m - order of the associated Legendre functions

   result - a double pointer to an array of (2*bw) samples
            representing the evaluation of the Legendre
	    series at (2*bw) Chebyshev nodes.

   trans_cos_pml_table - double pointer to array representing
                         the linearized form of trans(P) above.
			 See cospmls.{h,c} for a description
			 of the function Transpose_CosPmlTableGen()
			 which generates this array.

   sin_values - when m is odd, need to factor in the sin(x) that
                is factored out of the generation of the values
		in trans(P).

   workspace - a double array of size 2*bw -> temp space involving
          intermediate array

   fplan - pointer to fftw plan with input array being fcos
           and output being result; I'll probably be using the
	   guru interface to execute - that way I can apply the
	   same plan to different arrays; the plan should be
   
           fftw_plan_r2r_1d( 2*bw, fcos, result,
	        	     FFTW_REDFT01, FFTW_ESTIMATE );

*/

void InvSemiNaiveReduced(double *coeffs,
			 int bw, 
			 int m, 
			 double *result, 
			 double *trans_cos_pml_table, 
			 double *sin_values,
			 double *workspace,
			 fftw_plan *fplan )
{
  double *trans_tableptr;
  double *assoc_offset;
  int i, j, rowsize;
  double *p;
  double *fcos, fcos0, fcos1, fcos2, fcos3;
  double fudge ;

  fcos = workspace ;

  /* for paranoia, zero out arrays */
  memset( fcos, 0, sizeof(double) * 2 * bw );
  memset( result, 0, sizeof(double) * 2 * bw );

  trans_tableptr = trans_cos_pml_table;
  p = trans_cos_pml_table;

  /* main loop - compute each value of fcos

  Note that all zeroes have been stripped out of the
  trans_cos_pml_table, so indexing is somewhat complicated.
  */

  for (i=0; i<bw; i++)
    {
      if (i == (bw-1))
	{
	  if ( m % 2 )
	    {
	      fcos[bw-1] = 0.0;
	      break;
	    }
	}

      rowsize = Transpose_RowSize(i, m, bw);
      if (i > m)
	assoc_offset = coeffs + (i - m) + (m % 2);
      else
	assoc_offset = coeffs + (i % 2);

      fcos0 = 0.0 ; fcos1 = 0.0; fcos2 = 0.0; fcos3 = 0.0;
	  
      for (j = 0; j < rowsize % 4; ++j)
	fcos0 += assoc_offset[2*j] * trans_tableptr[j];
	  
      for ( ; j < rowsize; j += 4){
	fcos0 += assoc_offset[2*j] * trans_tableptr[j];
	fcos1 += assoc_offset[2*(j+1)] * trans_tableptr[j+1];
	fcos2 += assoc_offset[2*(j+2)] * trans_tableptr[j+2];
	fcos3 += assoc_offset[2*(j+3)] * trans_tableptr[j+3];
      }
      fcos[i] = fcos0 + fcos1 + fcos2 + fcos3 ;

      trans_tableptr += rowsize;
    }
    

  /*
    now we have the cosine series for the result,
    so now evaluate the cosine series at 2*bw Chebyshev nodes 
  */

  /* scale coefficients prior to taking inverse DCT */
  fudge = 0.5 / sqrt((double) bw) ;
  for ( j = 1 ; j < 2*bw ; j ++ )
    fcos[j] *= fudge ;
  fcos[0] /= sqrt(2. * ((double) bw));

  /* now take the inverse dct */
  /* NOTE that I am using the guru interface */
  fftw_execute_r2r( *fplan,
		    fcos, result );

  /* if m is odd, then need to multiply by sin(x) at Chebyshev nodes */
  if ( m % 2 )
    {
      for (j=0; j<(2*bw); j++)
	result[j] *= sin_values[j];
    }

  trans_tableptr = p;

  /* amscray */

}

/************************************************************************/

/* SemiNaiveReduced computes the Legendre transform of data.
   This function has been designed to be a component in
   a full spherical harmonic transform.  

   data - pointer to double array of size (2*bw) containing
          function to be transformed.  Assumes sampling at Chebyshev nodes

   bw   - bandwidth of the problem
   m   - order of the problem.  0 <= m < bw

   result - pointer to double array of length bw for returning computed
            Legendre coefficients.  Contains 
            bw-m coeffs, with the <f,P(m,m)> coefficient
            located in result[0]

   cos_pml_table - a pointer to an array containing the cosine
                   series coefficients of the Pmls (or Gmls)
		   for this problem.  This table can be computed
		   using the CosPmlTableGen() function, and
		   the offset for a particular Pml can be had
		   by calling the function NewTableOffset().
		   The size of the table is computed using
		   the TableSize() function.  Note that
		   since the cosine series are always zero-striped,
		   the zeroes have been removed.

   weights -> ptr to double array of size 4*bw - this array holds
           the weights for both even (starting at weights[0])
	   and odd (weights[2*bw]) transforms


   workspace -> tmp space: ptr to double array of size 4*bw

   fplan -> pointer to fftw plan with input array being weighted_data
           and output being cos_data; I'll probably be using the
	   guru interface to execute; the plan should be

	   fftw_plan_r2r_1d( 2*bw, weighted_data, cos_data,
			FFTW_REDFT10, FFTW_ESTIMATE ) ;


*/

void SemiNaiveReduced(double *data, 
		      int bw, 
		      int m, 
		      double *result,
		      double *workspace,
		      double *cos_pml_table, 
		      double *weights,
		      fftw_plan *fplan )
{
  int i, j, n;
  double result0, result1, result2, result3;
  double fudge ;
  double d_bw;
  int toggle ;
  double *pml_ptr, *weighted_data, *cos_data ;

  n = 2*bw;
  d_bw = (double) bw;

  weighted_data = workspace ;
  cos_data = weighted_data + (2*bw) ;

  /* for paranoia, zero out the result array */
  memset( result, 0, sizeof(double)*(bw-m));
 
  /*
    need to apply quadrature weights to the data and compute
    the cosine transform
  */
  if ( m % 2 )
    for ( i = 0; i < n    ; ++i )
      weighted_data[i] = data[ i ] * weights[ 2*bw + i ];
  else
    for ( i = 0; i < n    ; ++i )
      weighted_data[i] = data[ i ] * weights[ i ];

  /*
    smooth the weighted signal
  */

  fftw_execute_r2r( *fplan,
		    weighted_data,
		    cos_data );

  /* need to normalize */
  cos_data[0] *= 0.707106781186547 ;
  fudge = 1./sqrt(2. * ((double) n ) );
  for ( j = 0 ; j < n ; j ++ )
    cos_data[j] *= fudge ;

  /*
    do the projections; Note that the cos_pml_table has
    had all the zeroes stripped out so the indexing is
    complicated somewhat
  */
  

  /******** this is the original loop

  toggle = 0 ;
  for (i=m; i<bw; i++)
  {
  pml_ptr = cos_pml_table + NewTableOffset(m,i);

  if ((m % 2) == 0)
  {
  for (j=0; j<(i/2)+1; j++)
  result[i-m] += cos_data[(2*j)+toggle] * pml_ptr[j];
  }
  else
  {
  if (((i-m) % 2) == 0)
  {
  for (j=0; j<(i/2)+1; j++)
  result[i-m] += cos_data[(2*j)+toggle] * pml_ptr[j];
  }
  else
  {
  for (j=0; j<(i/2); j++)
  result[i-m] += cos_data[(2*j)+toggle] * pml_ptr[j];
  }
  } 
      
  toggle = (toggle+1) % 2;
  }

  *****/
 
  /******** this is the new loop *********/
  toggle = 0 ;
  for ( i=m; i<bw; i++ )
    {
      pml_ptr = cos_pml_table + NewTableOffset(m,i);

      result0 = 0.0 ; result1 = 0.0 ;
      result2 = 0.0 ; result3 = 0.0 ; 

      for ( j = 0 ; j < ( (i/2) % 4 ) ; ++j )
	result0 += cos_data[(2*j)+toggle] * pml_ptr[j];

      for ( ; j < (i/2) ; j += 4 )
	{
	  result0 += cos_data[(2*j)+toggle] * pml_ptr[j];
	  result1 += cos_data[(2*(j+1))+toggle] * pml_ptr[j+1];
	  result2 += cos_data[(2*(j+2))+toggle] * pml_ptr[j+2];
	  result3 += cos_data[(2*(j+3))+toggle] * pml_ptr[j+3];
	}

      if ((((i-m) % 2) == 0 ) || ( (m % 2) == 0 ))
	result0 += cos_data[(2*(i/2))+toggle] * pml_ptr[(i/2)];

      result[i-m] = result0 + result1 + result2 + result3 ;
	  
      toggle = (toggle + 1)%2 ;
	  
    }
}


