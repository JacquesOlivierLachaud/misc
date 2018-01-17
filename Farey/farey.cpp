/**
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 **/

/**
 * @file farey.cpp
 * @ingroup Examples
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5807), University of Savoie, France
 * @date 2015/10/30
 *
 * An example file named farey.cpp
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/arithmetic/LightSternBrocot.h"
#include "DGtal/arithmetic/LighterSternBrocot.h"

///////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace DGtal;
using namespace DGtal::Z2i;
///////////////////////////////////////////////////////////////////////////////
template <typename Fraction>
struct Farey
{
  typedef typename Fraction::Integer Integer;
  typedef std::list<Fraction> Sequence;
  typedef typename std::list<Fraction>::const_iterator ConstIterator;
  typedef typename std::list<Fraction>::iterator Iterator;

  Integer myN;
  Sequence mySequence;
  std::vector<Integer> myQuotientNbs;
  std::vector<Integer> myQuotientSums;

  Farey() 
    : myN( 1 )
  {
    addFraction( mySequence.end(), Fraction( Integer( 0 ), Integer( 1 ) ) );
    addFraction( mySequence.end(), Fraction( Integer( 1 ), Integer( 1 ) ) );
  }

  void addFraction( Iterator it, Fraction f )
  {
    mySequence.insert( it, f );
    for ( typename Fraction::ConstIterator itF = f.begin(), itFEnd = f.end();
          itF != itFEnd; ++itF )
      {
        Integer k   = (*itF).second;
        Integer u_k = (*itF).first;
        if ( k >= myQuotientNbs.size() ) 
          {
            myQuotientNbs.push_back( Integer( 1 ) );
            myQuotientSums.push_back( u_k );
          }
        else
          {
            myQuotientNbs[ k ]  += Integer( 1 );
            myQuotientSums[ k ] += u_k;
          }
      }
  }

  Integer       n() const     { return myN; }
  Integer       size() const  { return mySequence.size(); };
  ConstIterator begin() const { return mySequence.begin(); }
  ConstIterator end() const   { return mySequence.end(); }
  Integer       next() 
  {
    Integer nb = 0;
    myN       += 1;
    Iterator it = mySequence.begin(); ++it;
    Iterator itE = mySequence.end();  --itE;
    if ( it == itE ) 
      {
        addFraction( it, Fraction( Integer( 1 ), Integer( 2 ) ) );
        nb += 1;
      }
    while ( it != itE )
      {
        Iterator itNext = it; ++itNext;
        Fraction nleft = it->left();
        // std::cout << " " << nleft.p() << "/" << nleft.q() << std::endl;
        if ( nleft.q() == myN ) 
          {
            addFraction( it, nleft );
            nb += 1;
          }
        Fraction nright = it->right();
        if ( nright.q() == myN )
          {
            addFraction( itNext, nright );
            nb += 1;
          }
        it = itNext;
      }
    return nb;
  }
  
  Integer getSumOfQuotients() const
  {
    Integer s = 0;
    for ( typename std::vector<Integer>::const_iterator it = myQuotientSums.begin(), 
            itE = myQuotientSums.end(); it != itE; ++it )
      s += *it;
    return s;
  }

  Integer getNbOfQuotients() const
  {
    Integer s = 0;
    for ( typename std::vector<Integer>::const_iterator it = myQuotientNbs.begin(), 
            itE = myQuotientNbs.end(); it != itE; ++it )
      s += *it;
    return s;
  }
  int depth() const 
  {
    return myQuotientNbs.size() - 1;
  }
};

int main( int argc, char* argv[] )
{
  typedef DGtal::int64_t Integer;
  typedef LighterSternBrocot<Integer,Integer,DGtal::StdMapRebinder>::Fraction Fraction; // arbitrary large fractions

  Integer n = argc > 1 ? atoi( argv[ 1 ] ) : 100;

  Farey<Fraction> F;
  while ( F.n() != n )
    {
      std::cout << "- # F(" << F.n() << ") = " << F.size() 
                << " S_u = " << F.getSumOfQuotients()
                << " N_u = " << F.getNbOfQuotients()
                << " Avg_u = " << (double) F.getSumOfQuotients() / (double) F.getNbOfQuotients()
                << std::endl;
      std::cout << "  Nb_u[] =";
      for ( int i = 0; i <= F.depth(); ++i )
        {
          std::cout << " " << F.myQuotientNbs[ i ];
        }
      std::cout << endl;
      std::cout << " Avg_u[] =";
      for ( int i = 0; i <= F.depth(); ++i )
        {
          std::cout << " " << (double) F.myQuotientSums[ i ] / (double) F.myQuotientNbs[ i ];
        }
      std::cout << endl;
      // for ( Farey<Fraction>::ConstIterator it = F.begin(), itE = F.end();
      //       it != itE; ++it )
      //   std::cout << " " << it->p() << "/" << it->q();
      // std::cout << std::endl;
      Integer phi = F.next();
      std::cout << "  phi(" << F.n() << ") = " << phi << std::endl;
    }
  return 0;
}
