%% Generate a mixed grid(includeing triangles and quadrilaterals) on a part of a sphere
%% phi:   pi/4 ---> 3*pi/4
%% theta: 0    ---> pi/2

N = 2;
h = pi/2/N;


nodes = zeros((N+1)*(N+1), 3);
phi = pi/4 + [0:h:pi/2]; 
theta = [0:h:pi/2];

map_node = @(i,j) 

for i = 1 : N/2
	

	


end









