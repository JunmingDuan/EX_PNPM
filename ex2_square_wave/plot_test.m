function plot_ex1(K, n);

format long;
hold on;

%numer1 = load(['ex1_Nx',num2str(n),'_K',num2str(K),'.dat']);
%x1 = numer1(:,1); y1 = numer1(:,3); y2 = numer1(:,4);
%plot(x1, y1, 'o', x1, y2, '-k');
%legend(['k=',num2str(K)], 'exact');

for K = 2:4;
  numer1 = load(['ex1_Nx',num2str(n),'_K',num2str(K),'.dat']);
  x1 = numer1(:,1); y1 = numer1(:,3); y2 = numer1(:,4);
  plot(x1, y1, '-');
end
plot(x1, y2, '-k');
for K = 2:4;
  numer1 = load(['ex1_Nx',num2str(n),'_K',num2str(K),'_wrong.dat']);
  x1 = numer1(:,1); y1 = numer1(:,3); y2 = numer1(:,4);
  plot(x1, y1, '-');
end
legend('k=1', 'k=2', 'k=3', 'exact', 'k=1 wrong', 'k=2 wrong', 'k=3 wrong');

