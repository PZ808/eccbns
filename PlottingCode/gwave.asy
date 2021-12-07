include graph;
 
real width=15cm; 
real aspect=0.4; 
 
picture pic0, pic1; 

size(pic0, width, aspect*width, IgnoreAspect); 
size(pic1, width, aspect*width, IgnoreAspect); 
 
scale(pic0, true); 
scale(pic1, true); 
 
real xmin0 = 0; 
real xmax0 = 10;
real xmin1 = 0;
real xmax1 = 4; 

xlimits(pic0, 0, 6);
xlimits(pic1, 0, 6);

ylimits(pic0, -3, 3);
ylimits(pic1, -3, 3);
 
/* load the GW data from the files */ 
file data0 = line(input("xwave_3_4_1.4_10.0_0.00_1.0e-16.txt"));
file data1 = line(input("xwave_3_4_1.4_10.0_0.40_1.0e-12.txt"));

real[][] a = dimension(data0,0,0);
real[][] b = dimension(data1,0,0);

a = transpose(a);
b = transpose(b);

real[] t0 = a[1];
real[] h0 = a[4];

real[] t1 = b[1];
real[] h1 = b[4];

path circWave = graph(t0,h0);
path eccWave  = graph(t1,h1);

draw(pic0, circWave, black);
draw(pic1, eccWave, black);

xaxis(pic0, "", BottomTop, LeftTicks);
yaxis(pic0, "\LARGE{$h(t)$}", LeftRight, RightTicks);

xaxis(pic1, "\LARGE{$t$}", BottomTop, LeftTicks);
yaxis(pic1, "\LARGE{$h(t)$}", LeftRight, RightTicks);

yequals(pic0,0,Dotted); 
yequals(pic1,0,Dotted); 
 
pair min0=point(pic0,SW); 
pair max0=point(pic0,NE); 

pair min1=point(pic1,SW); 
pair max1=point(pic1,NE); 

real scale = (max0.x - min0.x) / (max1.x - min1.x); 
real shift = min0.x/scale - min0.x; 
 
transform t0 = pic0.calculateTransform(); 
transform t1 = pic1.calculateTransform(); 
transform T  = xscale(scale*t0.xx) * yscale(t1.yy); 
 
add(pic0.fit()); 

label("\LARGE{$e=0$}", 80NW+100E, UnFill(1mm));

real height = truepoint(N).y - truepoint(S).y; 

add(shift(0,-height)*(shift(shift)*pic1).fit(T)); 

//label("\large{$e=0.4$}", truepoint(SW) + 350E + 150N, UnFill(1mm));
label("\large{$e=0.4$}", truepoint(SW) + 80E + 160N, UnFill(1mm));

