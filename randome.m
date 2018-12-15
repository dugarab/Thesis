clc;
clear all;
close all;

Area = 1; %area of half circle
N = 1000;
pts_in = 0;
pts_out = 0;
x = zeros(1,N);
y = zeros(1,N);
x_in = zeros(1,N);
y_in = zeros(1,N);
x_out = zeros(1,N);
y_out = zeros(1,N);
for i= 1:N

x(i) = 6*rand-3; % [-3,3]
y(i) = 6*rand-3; %[-3,3]

if 3 >= sqrt(x(i)^2+ y(i)^2)
pts_in = pts_in +1;
x_in(pts_in) = x(i);
y_in(pts_in) = y(i);

else

pts_out = pts_out + 1;
x_out(pts_out) = x(i);
y_out(pts_out) = y(i);
end
end

errorpct= pts_out/1000;
disp(errorpct)
plot(x_in,y_in,'o', x_out,y_out,'+');
Aprprox_int = Area * (pts_in / N)