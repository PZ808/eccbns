include graph;
 
real width=12cm; 
real height=12cm;
real aspect=0.4; 

picture trajPic;

size(trajPic, width, width,IgnoreAspect);

// set the scale
scale(trajPic, true);

xlimits(trajPic, -40, 40);
ylimits(trajPic, -40, 40);


// get the data from the lines in the file           
file traject_data = line(input("xtraj_3_4_1.4_10.0_0.40_1.0e-12.txt"));

// read in the data into an array
real[][] coords = dimension(traject_data, 0, 0);

// transpose the array
coords = transpose(coords);

// set the coordinates 
real[] x = coords[1];
real[] y = coords[2];

path xy_traj = graph(x, y);

draw(trajPic, xy_traj, blue);

// axes 
xaxis(trajPic, "\large{$x\,(M_\odot)$}", BottomTop, LeftTicks);
yaxis(trajPic, "\large{$y\,(M_\odot)$}", LeftRight, RightTicks);  

add(trajPic.fit());
