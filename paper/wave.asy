include graph;
 
real width=15cm; 
real aspect=0.3; 
 
picture pic0, pic1,pic2;  

size(pic0,width,aspect*width,IgnoreAspect); 
size(pic1,width,aspect*width,IgnoreAspect); 
size(pic2,width,aspect*width,IgnoreAspect); 
 
scale(pic0,false); 
scale(pic1,false); 
scale(pic2,false); 
 
real xmin0=0; 
real xmax0=10;
real xmin1=0; 
real xmax1=10;
real xmin2=0;
real xmax2=3; 
 
/* load the GW data from the files */ 
file data0 = line(input("xwave_3_4_1.4_10.0_0.00_1.0e-16.txt"));
file data1 = line(input("xwave_3_4_1.4_10.0_0.10_1.0e-08.txt"));
file data2 = line(input("xwave_3_4_1.4_10.0_0.40_1.0e-12.txt"));

real[][] a = dimension(data0,0,0);
real[][] b = dimension(data1,0,0);
real[][] c = dimension(data2,0,0);

a = transpose(a);
b = transpose(b);
c = transpose(c);

real[] t0 = a[1];
real[] h0 = a[4];

real[] t1 = b[1];
real[] h1 = b[4];

real[] t2 = c[1];
real[] h2 = c[4];

path circWave = graph(t0,h0);
path eccWave1 = graph(t1,h1);
path eccWave2 = graph(t2,h2);

draw(pic0, circWave, red);
draw(pic1, eccWave1, blue);
draw(pic2, eccWave2, green);

//xaxis(pic0, "$t$", BottomTop, LeftTicks);
yaxis(pic0, "$h(t)$", LeftRight, RightTicks);
xaxis(pic1, "$t$", BottomTop, LeftTicks);
yaxis(pic1, "$h(t)$", LeftRight, RightTicks);
xaxis(pic2, "$t$", BottomTop, LeftTicks);
yaxis(pic2, "$h(t)$", LeftRight, RightTicks);

yequals(pic0,0,Dotted); 
yequals(pic1,0,Dotted); 
yequals(pic2,0,Dotted); 
 
pair min0=point(pic0,SW); 
pair max0=point(pic0,NE); 
pair min1=point(pic1,SW); 
pair max1=point(pic1,NE); 
pair min2=point(pic2,SW); 
pair max2=point(pic2,NE); 

real scale=(max1.x-min1.x)/(max2.x-min2.x); 
real shift=min1.x/scale-min2.x; 
 
transform t1 = pic1.calculateTransform(); 
transform t2 = pic2.calculateTransform(); 
transform T=xscale(scale*t1.xx)*yscale(t2.yy); 
 
add(pic0.fit()); 
//add(pic1.fit()); 

real height=truepoint(N).y-truepoint(S).y; 

//add(shift(0,-height)*(shift(shift)*pic1).fit(T)); 
add(shift(0,-height)*(shift(shift)*pic2).fit(T)); 

