function plot_ex1(K, n);

format long;
hold on;
mu = 1;
plot(linspace(0,1,1000), sin(2*pi*(linspace(0,1,1000)-mu*3e-3)), '-k');

numer1 = load(['ex1_Nx',num2str(n),'_K',num2str(K),'.dat']);
x1 = numer1(:,1); y1 = numer1(:,3);
plot(x1, y1, 'o');

ind = 18;
plot(x1(K*ind+1:K*(ind+1)), y1(K*ind+1:K*(ind+1)), '*');

