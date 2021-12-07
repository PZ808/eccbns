import graph;
import interpolate;

// We want to interpolate the data 
picture pic0, pic1;


size(pic0,width,aspect*width,IgnoreAspect); 
size(pic1,width,aspect*width,IgnoreAspect);  


//size(135mm,100mm, IgnoreAspect);

//scale(false, true);

string fname0="";
string fname1="";
string fname2="";

pen p=linewidth(1);
file in0=line(input(fname0);
file in1=line(input(fname1);
//file in2=line(input(fname2);

real[][] a=dimension(in0,0,0);
real[][] b=dimension(in1,0,0);
//real[][] c=dimension(in2,0,0);

a=transpose(a);
b=transpose(b);
//c=transpose(c);

real[] e0 = a[0];
real[] e1 = b[0];
//real[] e2 = c[0];

real[] olap0 = a[1];
real[] olap1 = b[1];
//real[] olap2 = c[1];

draw( 
      pic0,
      graph( e0, olap1, Hermite(natural) ), 
      red, 
      "$M = 2.8M_{\odot}$" 
    ); 

xaxis( "$e_0$",BottomTop, LeftTicks );
yaxis( "overlap",LeftRight,RightTicks );             

file in2=line( input("") );
real[][] b=dimension( in2, 0 ,0 );
b=transpose(b);
real[] e2 = b[0];
real[] olap2 = b[1];

draw( 
      graph(e2,olap2), 
      dashed+black+p, 
      "$M = 6.4M_{\odot}$"
    ); 

//ylimits(0.9963,1.001);

attach( legend(), truepoint(E),10E, UnFill );
