% x = β and y = γ
f = @(x,y) -(x/N)*S(t)*I(t);
f = @(x,y) (x/N)*S(t)*I(t)-y*I(t);
f = @(x,y) y*I(t);
h = 1;
t = (0:h:100);
S(0) = 990;
I(0) = 10;
R(0) = 0;
N = S(t)+I(t)+R(t);
x=0.3;
y=0.1;
for i = 1:(length(t)-1)
k1 = f(x(i),y(i));
k2 = f(x(i)+.5*h,y(i)+.5*k1*h);
k3 = f(x(i)+.5*h,y(i)+.5*k2*h);
k4 = f(x(i)+h,y(i)+k3*h);
y(i+1) = y(i)+(1/6)*(k1+2*k2+2*k3+k4)*h;

plot (x,y);
xlabel('β');
ylabel('γ');
title('4th-Order Runge-Kutta Method');
end