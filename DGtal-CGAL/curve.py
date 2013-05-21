import math;

# helicoid
# def f(t):
#     return (math.cos(5*t),math.sin(5*t),math.cos(t))
def f(t):
    return (0.2*t*math.cos(5*t),(1-0.2*t)*math.sin(5*t),0.3*t*math.cos(t))

def output(h):
    b = int( 5.0*2.0*3.14159626*math.ceil(1.0/h) );
    for i in range(0,b):
        t = ( h*float( i )/5.0 )
        (x,y,z)=f( t )
        print int( math.floor(x/h) ), int( math.floor(y/h) ), int( math.floor(z/h) )

output(0.02)

