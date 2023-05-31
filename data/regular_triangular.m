clear;
xmin = 0; xmax = 1; ymin = 0; ymax = 1;
N = 11;
x = linspace(xmin, xmax, N);
y = linspace(ymin, ymax, N);
p = zeros(N^2, 2);
for i = 1:N
	p((i-1)*N+1:i*N, 1) = x(i);
	p((i-1)*N+1:i*N, 2) = y;
end
%plot(p(:,1), p(:,2), 'x')

dt = delaunayTriangulation(p);
edge_list = edges(dt);

%triplot(dt, 'color', 'black'); 
%axis square
% export data
writematrix(dt.Points, 'point_list.txt', 'Delimiter', '\t');
writematrix(dt.ConnectivityList, 'connectivity_list.txt', 'Delimiter', '\t')
writematrix(edge_list, 'edge_list.txt', 'Delimiter', '\t')