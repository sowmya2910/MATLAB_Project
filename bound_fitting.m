function bound_fitting
th1_lw = -pi/2;
th1_up = pi/3;

th2_lw = -10*pi/180;
th2_up = 10*pi/180;

th3_lw = -pi/3;
th3_up = pi/3;

x = -3*pi:pi/256:3*pi;
cost = 1E5;
y1 = cost*(x<th1_lw | x>th1_up);
y2 = cost*(x<th2_lw | x>th2_up);
y3 = cost*(x<th3_lw | x>th3_up);
y1f = polyfit(x,y1,14);
y2f = polyfit(x,y2,14);
y3f = polyfit(x,y3,14);
save('costfilter.mat','y1f','y2f','y3f')