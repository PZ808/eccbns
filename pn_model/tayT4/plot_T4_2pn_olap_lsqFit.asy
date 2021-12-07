import graph;
import interpolate;

size(12cm, 10cm, IgnoreAspect);

pen p=linewidth(1.1);

marker mark=marker( scale(1mm)*cross(4, false, r=0.35), 
                    blue, Fill ); 

ylimits( 0.55, 1.00 );

file in1=line(input("Match_T4vsXmodel_4_1.4_1.4_3c4rPN.txt"));
real[][] a=dimension(in1,0,0);
a=transpose(a);
real[] e=a[0];
real[] olap=a[1];

draw( 
      graph(e, olap, Hermite(natural)), 
      red + p, 
      "$M = 2.8\,M_{\odot}$" 
    ); 

xaxis("\LARGE{$e_0$}",BottomTop, LeftTicks);
yaxis("\large{overlap}",LeftRight,RightTicks);             

file in2=line(input("Match_T4vsXmodel_4_1.4_5.0_3c4rPN.txt"));
real[][] b=dimension(in2, 0 ,0);
b=transpose(b);
real[] e2 = b[0];
real[] olap2 = b[1];

draw( 
      graph(e2, olap2, Hermite(natural)), 
      dashed + blue + p, 
      "$M = 6.4\,M_{\odot}$"
    ); 

file in3=line(input("Match_T4vsXmodel_4_5.0_5.0_3c4rPN.txt"));
real[][] c=dimension(in3, 0 ,0);
c=transpose(c);
real[] e3 = c[0];
real[] olap3 = c[1];

draw( 
      graph(e3, olap3, Hermite(natural)), 
      linetype("8 8 0 8") + p,
      "$M = 10.0\,M_{\odot}$"
      //mark 
    ); 

attach( legend(), point(SW), 50N + 35E, UnFill );
