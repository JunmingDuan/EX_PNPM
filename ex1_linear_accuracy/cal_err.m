function plot_ex1(K, n);

format long;
hold on;
mu = 1;
t = 0.5;

numer1 = load(['ex1_Nx',num2str(n),'_K',num2str(K),'.dat']);
x1 = numer1(:,1); y1 = numer1(:,3); y2 = numer1(:,4);
plot(x1, y1, 'o-', x1, y2, '-r');

m = 6;
err = zeros(m, 3);
for i = 1:m;
  n = 10*2^(i-1);
  %n = 2^(i+1);
  numer1 = load(['ex1_Nx',num2str(n),'_K',num2str(K),'.dat']);
  w1 = numer1(:,2); y1 = numer1(:,3); y2 = numer1(:,4);
  err(i, 1) = sqrt(sum((y1-y2).^2.*w1)/n);
  err(i, 2) = sum(abs(y1-y2).*w1)/n;
  err(i, 3) = max(abs(y1-y2));
end
order = zeros(m, 3);
for i = 1:m-1;
  order(i+1,:) = log2(err(i,:)./err(i+1,:));
end
%err
%order

diary table1.dat
diary on;
for j = 1:m;
  fprintf('%3d ', 10*2^(j-1));
  for i = 1:3
    fprintf('& %.3e & %.2f ', err(j,i), order(j,i));
  end
  fprintf('\\\\ \n');
  %fprintf('%.3e\n', min_y(j));
end
diary off;


