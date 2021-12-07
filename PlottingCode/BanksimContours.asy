import graph;
import contour;
import palette;

size(14cm, 14cm, IgnoreAspect);

defaultpen(1bp);
int Divs = 10;
real Low = 0.4;
real High = 1.0;                   

xlimits(0,0.4);

pen Tickpen = black;
pen[] Palette = Rainbow();

scale(false);

file in = line(input("bankRes35PN_10025INJecc_mt_fix02.dat"));
real[][] data = dimension(in, 0, 0);
data = transpose(data);

real[] ecc  = data[0];
real[] mass = data[1];
real[] olap = data[2];

int n=ecc.length;

pair[] points = new pair[n];
real[] vals   = new real[n];

/* major contours */
real[] Cvals=uniform(Low, High, Divs);

for (int i=0; i<n; ++i) {
  points[i] = (ecc[i], mass[i]);
  vals[i] = olap[i];
}     

draw( contour(points, vals, new real[]{0.7, 0.8, 0.9, 0.95}, operator --), Tickpen );
//draw(contour(points, vals, Cvals, operator --), Tickpen );

xaxis("$e_0$", BottomTop, LeftTicks, above=true);
yaxis("$M$", LeftRight, RightTicks, above=true);

picture bar;

bounds range = image(points, vals, Palette);

palette(bar,"$\bar\mathcal{F}$", range,(0,0), (0.9cm,12cm), Right, Palette,
        PaletteTicks("$%+#.1f$"));

add(bar.fit(),point(E),30E);

