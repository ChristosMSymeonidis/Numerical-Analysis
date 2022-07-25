clear;
clf;
close all;
syms x
m=2*4;
k=3*4;
j=1;
c=1;
f = @(x) (cos(m*x)+sin(k*x)).^2;
int = integral(f,-1,1);

for N = [100 500 1000 2000] 
count = 0;
for i = 1:N    
xr=rand*2-1;
yr=rand*4;
if yr<=f(xr)    
    inside(i,1) = xr;
    inside(i,2) = yr;
    count = count + 1;
        else    
    outside(i,1) = xr;
    outside(i,2) = yr;
end
end
MCI = count*8/N;

subplot(2,4,j);
plot (inside(:,1),inside(:,2),'.b');
hold on
plot (outside(:,1),outside(:,2),'.r');
hold on
fplot (f(x),[-1,1],'k');
title("MCI ("+N+") = "+MCI+" ("+int+")");
grid on 
j=j+1;
end


for N=10:10:2000
    tic;
    count = 0;
for i = 1:N    
xr=rand*2-1;
yr=rand*4;
if yr<=f(xr)    
    inside(i,1) = xr;
    inside(i,2) = yr;
    count = count + 1;
        else    
    outside(i,1) = xr;
    outside(i,2) = yr;
end
end
MCi = count*8/N;
time(c) = toc;   
error(c) = abs(1-(MCi/int))*100;
c=c+1;
end

subplot(2,4,5:6);
plot(10:10:2000, time,'- .b');
title("Time - Number of Samples");
ylabel("Time"); xlabel("Number of Samples");
grid on;

subplot(2,4,7:8);
plot(10:10:2000, error,'- .r');
title("%Error - Number of Samples");
ylabel("%Error"); xlabel("Number of Samples")
grid on;