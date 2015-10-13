// gcc -std=c99 -Wall -pedantic hierarchy.c -lm -o hierarchy
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

// Not in math.h
#define PI 3.14159265358979323846
#define SIZE 256
#define LVL 8

typedef float Value;

Value square_int( int x );
Value square( Value x );
Value distance2( int x1, int y1, int z1, int x2, int y2, int z2 );
Value distance( int x1, int y1, int z1, int x2, int y2, int z2 );

typedef Value (* VoxelFunctor )( Value data, int x, int y, int z );
Value moment000( Value data, int x, int y, int z );
Value moment100( Value data, int x, int y, int z );
Value moment010( Value data, int x, int y, int z );
Value moment001( Value data, int x, int y, int z );
Value moment200( Value data, int x, int y, int z );
Value moment020( Value data, int x, int y, int z );
Value moment002( Value data, int x, int y, int z );
Value moment110( Value data, int x, int y, int z );
Value moment101( Value data, int x, int y, int z );
Value moment011( Value data, int x, int y, int z );

struct SImage {
  int size;
  int lvl;
  Value* data;
};
typedef struct SImage Image;

void  Image_init( Image* img, int k );
Value Image_get( Image* img, int x, int y, int z );
void  Image_set( Image* img, int x, int y, int z, Value v );
void  Image_finish( Image* img );

struct SMipMap {
  int   max_lvl;
  Image hierarchy[ LVL+1 ];
};
typedef struct SMipMap MipMap;

void   MipMap_init_from_image_and_functor( MipMap* mipmap, Image* img, VoxelFunctor f );
Image* MipMap_get_image( MipMap* mipmap, int lvl );
void   MipMap_finish( MipMap* mipmap );

int nb_iteration_exact     = 0;
int nb_access_exact        = 0;
int nb_iteration_hierarchy = 0;
int nb_access_hierarchy    = 0;
 

Value square_int( int x )
{
  Value f = (Value) x;
  return f*f;
}
Value square( Value x )
{
  return x*x;
}

Value distance2( int x1, int y1, int z1, int x2, int y2, int z2 )
{
  return square_int(x1 - x2) + square_int(y1 - y2) + square_int(z1 - z2);
}

Value distance( int x1, int y1, int z1, int x2, int y2, int z2 )
{
  return sqrt( distance2( x1, y1, z1, x2, y2, z2 ) );
}

void Image_init( Image* img, int k )
{
  img->lvl = k;
  img->size = 1 << k;
  img->data = (Value*) malloc( (img->size)*(img->size)*(img->size)*sizeof( Value ) );
}

Value Image_get( Image* img, int x, int y, int z )
{
  int k = img->lvl;
  return img->data[ ( ( ( z << k ) + y ) << k ) + x ];
}

void Image_set( Image* img, int x, int y, int z, Value v )
{
  int k = img->lvl;
  img->data[ ( ( ( z << k ) + y ) << k ) + x ] = v;
}

void Image_finish( Image* img )
{
  img->lvl = 0;
  img->size = 0;
  free( img->data );
  img->data = 0;
}

void Image_ball( Image* img, int x0, int y0, int z0, float r )
{
  Value r2 = r*r;
  Value* data = img->data;
  for ( int z = 0; z < img->size; ++z )
    for ( int y = 0; y < img->size; ++y )
      for ( int x = 0; x < img->size; ++x )
        {
          if ( distance2( x0, y0, z0, x, y, z ) <= r2 )
            *data++ = (Value) 1;
          else
            *data++ = (Value) 0;
        }
}

void MipMap_init_from_image_and_functor( MipMap* mipmap, Image* img, VoxelFunctor f )
{
  assert( ( img->lvl >= 1 ) && ( img->lvl <= LVL ) );
  mipmap->max_lvl = img->lvl;
  for ( int k = img->lvl; k > 0; --k )
    {
      Image_init( &mipmap->hierarchy[ k ], k );
    }
  // Copy image with functor (a moment).
  Image* src = img;
  Image* dst = MipMap_get_image( mipmap, img->lvl );
  for ( int z = 0; z < src->size; ++z )
    for ( int y = 0; y < src->size; ++y )
      for ( int x = 0; x < src->size; ++x )
        Image_set( dst, x, y, z, f( Image_get( src, x, y, z ), x, y, z ) );
  // Compute hierarchy
  for ( int k = img->lvl - 1; k > 0; --k )
    {
      Image* src = MipMap_get_image( mipmap, k+1 );
      Image* dst = MipMap_get_image( mipmap, k );
      for ( int z = 0; z < dst->size; ++z )
        for ( int y = 0; y < dst->size; ++y )
          for ( int x = 0; x < dst->size; ++x )
            {
              Value f = Image_get( src,   2*x,   2*y,   2*z );
              f      += Image_get( src, 2*x+1,   2*y,   2*z );
              f      += Image_get( src,   2*x, 2*y+1,   2*z );
              f      += Image_get( src, 2*x+1, 2*y+1,   2*z );
              f      += Image_get( src,   2*x,   2*y, 2*z+1 );
              f      += Image_get( src, 2*x+1,   2*y, 2*z+1 );
              f      += Image_get( src,   2*x, 2*y+1, 2*z+1 );
              f      += Image_get( src, 2*x+1, 2*y+1, 2*z+1 );
              Image_set( dst, x, y, z, f / (Value) 8 );
            }
    }
}

Image* MipMap_get_image( MipMap* mipmap, int lvl )
{
  assert( ( lvl >= 1 ) && ( lvl <= LVL ) );
  return &mipmap->hierarchy[ lvl ];
}

void   MipMap_finish( MipMap* mipmap )
{
  for ( int k = 1; k <= mipmap->max_lvl; ++k )
    Image_finish( MipMap_get_image( mipmap, k ) );
  mipmap->max_lvl = 0;
}

Value moment000( Value data, int x, int y, int z )
{
  return (Value) data;
}
Value moment100( Value data, int x, int y, int z )
{
  return (Value) data*x;
}
Value moment010( Value data, int x, int y, int z )
{
  return (Value) data*y;
}
Value moment001( Value data, int x, int y, int z )
{
  return (Value) data*z;
}
Value moment200( Value data, int x, int y, int z )
{
  return (Value) data*x*x;
}
Value moment020( Value data, int x, int y, int z )
{
  return (Value) data*y*y;
}
Value moment002( Value data, int x, int y, int z )
{
  return (Value) data*z*z;
}
Value moment110( Value data, int x, int y, int z )
{
  return (Value) data*x*y;
}
Value moment101( Value data, int x, int y, int z )
{
  return (Value) data*x*z;
}
Value moment011( Value data, int x, int y, int z )
{
  return (Value) data*y*z;
}

Value computeExact( MipMap* M, int x0, int y0, int z0, Value r )
{
  Image* img = MipMap_get_image( M, M->max_lvl );
  int minx = x0 - (int) ceil( r );
  int miny = y0 - (int) ceil( r );
  int minz = z0 - (int) ceil( r );
  int maxx = x0 + (int) ceil( r );
  int maxy = y0 + (int) ceil( r );
  int maxz = z0 + (int) ceil( r );
  minx = minx <= 0 ? 0 : minx;
  miny = miny <= 0 ? 0 : miny;
  minz = minz <= 0 ? 0 : minz;
  maxx = maxx >= img->size ? img->size-1 : maxx;
  maxy = maxy >= img->size ? img->size-1 : maxy;
  maxz = maxz >= img->size ? img->size-1 : maxz;
  Value r2  = r*r;
  Value acc = (Value) 0;
  for ( int z = minz; z <= maxz; ++z )
    for ( int y = miny; y <= maxy; ++y )
      for ( int x = minx; x <= maxx; ++x )
        {
          if ( distance2( x0, y0, z0, x, y, z ) <= r2 )
            {
              acc += Image_get( img, x, y, z );
              nb_access_exact += 1;
            }
          nb_iteration_exact += 1;
        }
  return acc;
}

void goUp( int xyzk[] )
{
  xyzk[ 0 ] >>= 1;
  xyzk[ 1 ] >>= 1;
  xyzk[ 2 ] >>= 1;
  xyzk[ 3 ]  -= 1;
}
void goDown( int xyzk[] )
{
  xyzk[ 0 ] <<= 1;
  xyzk[ 1 ] <<= 1;
  xyzk[ 2 ] <<= 1;
  xyzk[ 3 ]  += 1;
}

/* The binomial tree is seen as son / 8-brother tree. */
/* return 0 when finished. */
int goNext( int xyzk[] )
{
  while ( ( xyzk[ 0 ] & 0x1 ) && ( xyzk[ 1 ] & 0x1 ) && ( xyzk[ 2 ] & 0x1 ) )
    goUp( xyzk );
  if ( xyzk[ 3 ] == 0 ) { xyzk[ 3 ] = -1; return 0; }
  if ( ( xyzk[ 0 ] & 0x1 ) == 0 )  
    xyzk[ 0 ] += 1; 
  else {
    xyzk[ 0 ] &= ~0x1;
    if ( ( xyzk[ 1 ] & 0x1 ) == 0 )
      xyzk[ 1 ] += 1;
    else {
      xyzk[ 1 ] &= ~0x1;
      xyzk[ 2 ] += 1;
    }
  }
  return 1;
}

Value computeHierarchy( MipMap* M, int x0, int y0, int z0, Value r )
{
  Value radii [ LVL+1 ];
  Value weight[ LVL+1 ];
  Value diag  [ LVL+1 ];
  // int k           = 0;          // current level in the hierarchy
  int max_k       = M->max_lvl; // number of levels in the hierarchy
  // int h           = max_k;      // height in hierarchy (max_k - k )
  // int size        = 1 << max_k; // size of each cell
  radii [ max_k ] = r;
  weight[ max_k ] = (Value) 1;
  diag  [ max_k ] = (Value) (sqrt(3.0)/2.0);
  for ( int i = max_k - 1; i >= 0; --i )
    {
      radii [ i ] = radii [ i+1 ] / (Value) 2;
      weight[ i ] = weight[ i+1 ] * (Value) 8;
      diag  [ i ] = diag  [ i+1 ] * (Value) 2;
    }
  int xyzk[ 4 ] = { 0, 0, 0, 0 };
  Value acc = 0;
  Value r2  = r*r;
  do 
    {
      nb_iteration_hierarchy += 1;
      int k = xyzk[ 3 ];  // current level in the hierarchy
      int h = max_k - k;  // height in hierarchy (max_k - k )
      Value dK2    = distance2( ( 2*xyzk[ 0 ] + 1) << h, ( 2*xyzk[ 1 ] + 1) << h, ( 2*xyzk[ 2 ] + 1 ) << h, 
                                2*x0+1, 2*y0+1, 2*z0+1 );
      /* Value d2     = distance2( xyzk[ 0 ] << h, xyzk[ 1 ] << h, xyzk[ 2 ] << h,  */
      /*                           x0, y0, z0 ); */
      Value d2     = dK2 / (Value) 4;
      Value delta2 = square( diag[ k ] );
      Value upper2 = ( r2 >= delta2 ) ? r2 - 2.0*r*diag[ k ] + delta2 : -1.0; // (r-diag/2)^2
      Value lower2 = r2 + 2.0*r*diag[ k ] + delta2; // (r+diag/2)^2
      //printf( "[%d] %d %d %d | l2=%f d2=%f u2=%f\n", k, xyzk[ 0 ], xyzk[ 1 ], xyzk[ 2 ], lower2, d2, upper2 );
      // Takes care of finest cells.
      if ( h == 0 ) 
        {
          if ( d2 <= r2 ) 
            { // cell is completely inside
              // printf("[%d] %d %d %d\n", k, xyzk[ 0 ], xyzk[ 1 ], xyzk[ 2 ] );
              acc += Image_get( MipMap_get_image( M, k ), 
                                xyzk[ 0 ], xyzk[ 1 ], xyzk[ 2 ] ); 
              nb_access_hierarchy += 1;
              goNext( xyzk );
            }
          else // cell is completely outside
            goNext( xyzk );
        }
      else
        { // cell is completely inside
          if ( d2 <= upper2 ) 
            {
              // printf("[%d] %d %d %d\n", k, xyzk[ 0 ], xyzk[ 1 ], xyzk[ 2 ] );
              acc += weight[ k ] * Image_get( MipMap_get_image( M, k ), 
                                              xyzk[ 0 ], xyzk[ 1 ], xyzk[ 2 ] );
              nb_access_hierarchy += 1;
              goNext( xyzk );
            }
          else if ( d2 > lower2 )
            // cell is completely outside
            goNext( xyzk );                 
          else goDown( xyzk );
        }
    }
  while ( xyzk[ 3 ] > 0 );
  return acc;
}

int main( int argc, char* argv[] )
{
  int lvl = argc > 1 ? atoi( argv[ 1 ] ) : 6;
  Value R = argc > 2 ? atof( argv[ 2 ] ) : 15.0;
  Value r = argc > 3 ? atof( argv[ 3 ] ) : 5.0;

  Image I;
  Image_init( &I, lvl );

  printf("---- Creating Image of size (%d)^3 with ball of radius %f ----\n", I.size, R );
  int x0 = 1 << (lvl-1); 
  int y0 = 1 << (lvl-1); 
  int z0 = 1 << (lvl-1); 
  Image_ball( &I, x0, y0, z0 + (int) floor( R ), R );

  printf("---- Creating MipMap for volume ----\n" );
  MipMap vol;
  MipMap_init_from_image_and_functor( &vol, &I, moment000 );

  printf("---- Computing exact volume of shape cap ball(%f) ----\n", (double) r );
  Value exact = computeExact( &vol, x0, y0, z0, r );
  printf("     - expected value       ~ %f (should be below, the greater R the closest)\n", 4.0*PI/6.0*r*r*r );
  printf("     - exact discrete value = %f, iter/access %d/%d\n", exact, nb_iteration_exact, nb_access_exact );
  Value hier = computeHierarchy( &vol, x0, y0, z0, r );
  printf("     - hier. discrete value = %f, iter/access %d/%d\n", hier, nb_iteration_hierarchy, nb_access_hierarchy );
 
  printf("---- Freeing memory ----\n" );
  MipMap_finish( &vol );
  Image_finish( &I );
  return 0;
}
