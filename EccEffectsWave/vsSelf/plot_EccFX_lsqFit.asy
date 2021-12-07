import graph;
import interpolate;

size(12cm, 10cm, IgnoreAspect);

pen p=linewidth(1.1);

marker mark = marker(scale(1mm)*cross(4, false, r=0.35), blue, Fill); 

ylimits(0.10, 1.0005);

file in1=line(input("MatchSelf_1.4_1.4_3c4rPN.txt"));
real[][] a=dimension(in1,0,0);
a=transpose(a);
real[] e=a[0];
real[] olap=a[1];

draw( 
      graph(e, olap, Hermite(natural)), 
      red+p, 
      "$M = 2.8\,M_{\odot}$" 
    ); 

xaxis( "\LARGE{$e_0$}", BottomTop, LeftTicks );
yaxis( "\large{overlap}", LeftRight,RightTicks );             

file in2=line( input("MatchSelf_1.4_5.0_3c4rPN.txt") );
real[][] b=dimension( in2, 0 ,0 );
b=transpose(b);
real[] e2 = b[0];
real[] olap2 = b[1];

draw( 
      graph( e2, olap2, Hermite(natural) ), 
      dashed + blue + p, 
      "$M = 6.4\,M_{\odot}$"
    ); 

file in3=line( input("MatchSelf_5.0_5.0_3c4rPN.txt") );
real[][] c=dimension( in3, 0 ,0 );
c=transpose(c);
real[] e3 = c[0];
real[] olap3 = c[1];

draw( 
      graph( e3, olap3, Hermite(natural) ), 
      longdashdotted+p,
      "$M = 10\,M_{\odot}$"
    ); 

/*
file in4=line( input("MatchSelf_1.4_10.0_3c4rPN.txt") );
real[][] d=dimension( in3, 0 ,0 );
d=transpose(d);
real[] e4 = c[0];
real[] olap4 = c[1];
draw( graph( e4, olap4, Hermite(natural) ), dashdotted+blue, "$M = 11.4\,M_{\odot}$"); 
*/

attach(legend(), NE,-210E + 36N, UnFill);
