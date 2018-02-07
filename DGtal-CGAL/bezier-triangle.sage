var('d0','e0')
var('d1','e1')
var('d2','e2')
var('a','b','c')

M = matrix(
  [[1,0,0,0,0,0],
   [0,1,0,0,0,0],
   [0,0,1,0,0,0],
   [-d0-e0,d0,0,-e0,e0,e0],
   [0,-d1-e1,d1,e1,-e1,e1],
   [d2,0,-d2-e2,e2,e2,-e2] ] )

B = vector( [a,b,c,0,0,0] )
X = M.inverse() * B
# (a, b, c, -1/2*a*((d0 + e0)/e0 + (d2 - (d0 + e0)*e2/e0)/e2) + 1/2*b*(d0/e0 + (d1 - d0*e1/e0 + e1)/e1) - 1/2*c*(d1/e1 - (d2 + e2)/e2), -1/2*b*d0/e0 - 1/2*a*(d2 - (d0 + e0)*e2/e0)/e2 + 1/2*c*(d2 + e2)/e2, 1/2*a*(d0 + e0)/e0 + 1/2*b*(d1 - d0*e1/e0 + e1)/e1 - 1/2*c*d1/e1)
f(a,b,c,d0,e0,d1,e1,d2,e2)=(a, b, c, -1/2*a*((d0 + e0)/e0 + (d2 - (d0 + e0)*e2/e0)/e2) + 1/2*b*(d0/e0 + (d1\
 - d0*e1/e0 + e1)/e1) - 1/2*c*(d1/e1 - (d2 + e2)/e2), -1/2*b*d0/e0 - 1/2*a*(d2 - (\
d0 + e0)*e2/e0)/e2 + 1/2*c*(d2 + e2)/e2, 1/2*a*(d0 + e0)/e0 + 1/2*b*(d1 - d0*e1/e0\
 + e1)/e1 - 1/2*c*d1/e1)

# Triangle coin bas gauche avec gradients orth sur les trois arêtes
f(a,b,c,0,1,1/2,-1,-1,1)

b = vector( f(a,b,c,0,1,1/2,-1,-1,1) )
var('r','s','t')
R(x,y)=1-x-y
S(x,y)=x
T(x,y)=y
P(x,y)=R(x,y)^2*b[0] + S(x,y)^2 * b[1] + T(x,y)^2 * b[ 2 ] + R(x,y)*S(x,y) * b[3] + S(x,y)*T(x,y) * b[ 4 ] + T(x,y)*R(x,y) * b[5]

Q(a,b,c,x,y)=a*(x + y - 1)^2 - 1/4*(2*a + b + c)*(x + y - 1)*x + b*x^2 - 1/4*(2*a + b + c)*(x + y - 1)*y + a*x*y + c*y^2

Q(1,0,0,x,y)
(x + y - 1)^2 - 1/2*(x + y - 1)*x - 1/2*(x + y - 1)*y + x*y

