%% generate a mixed structure mesh on [0,1]^2
clear;
clc;

N = 2;
h = 1/N;

mapnode = @(i,j) (i-1)*(N+1)+j;
nodes = zeros((N+1)*(N+1), 3);
for i = 1 : N+1
    for j = 1 : N+1
        nodes(mapnode(i,j), 1) = (i-1)*h;
        nodes(mapnode(i,j), 2) = (j-1)*h;
        nodes(mapnode(i,j), 3) = 0;
    end
end

% 2*N*(N+1)+N*(N/2)

lines = zeros(2*N*(N+1)+N*(N/2), 3); % first:  node one
                                     % second: node two
                                     % third:  mark
                                     %          0    interior
                                     %          1    dirichlet
mapline = @(i,j) 2*N*(i-1) + (2*j-1);
line = 1;
for i = 1 : N
    for j = 1 : N
        lines(line, 1) = mapnode(i,j);
        lines(line, 2) = mapnode(i+1,j);
        if (j == 1)
            lines(line, 3) = 1;
        else
            lines(line, 3) = 0;
        end
        line++;

        lines(line, 1) = mapnode(i,j);
        lines(line, 2) = mapnode(i,j+1);
        if (i == 1)
            lines(line, 3) = 1;
        else
            lines(line, 3) = 0;
        end
        line++;
    end
end

for i = 1 : N
    lines(line, 1) = mapnode(i, N+1);
    lines(line, 2) = mapnode(i+1, N+1);
    lines(line, 3) = 1;
    line++;
end
for j = 1 : N
    lines(line, 1) = mapnode(N+1, j);
    lines(line, 2) = mapnode(N+1, j+1);
    lines(line, 3) = 1;
    line++;
end

for i = N/2+1 : N
    for j = 1 : N
        lines(line, 1) = mapnode(i, j+1);
        lines(line, 2) = mapnode(i+1, j);
        lines(line, 3) = 0;
        line++;
    end
end
line--; % number of lines


elems = zeros(3*N*N/2, 5); % first:  node one
                              % second: node two
                              % third:  node three
                              % fourth: node four
                              % fifth:  type
                              %          2    triangle
                              %          3    quadrilateral
elem = 1;
for i = 1 : N/2
    for j = 1 : N
        elems(elem, 1) = mapnode(i,j);
        elems(elem, 2) = mapnode(i+1,j);
        elems(elem, 3) = mapnode(i+1,j+1);
        elems(elem, 4) = mapnode(i,j+1);
        elems(elem, 5) = 3;
        elem++;
    end
end
for i = N/2+1 : N
    for j = 1 : N
        elems(elem, 1) = mapnode(i,j);
        elems(elem, 2) = mapnode(i+1,j);
        elems(elem, 3) = mapnode(i,j+1);
        elems(elem, 5) = 2;
        elem++;

        elems(elem, 1) = mapnode(i+1,j);
        elems(elem, 2) = mapnode(i+1,j+1);
        elems(elem, 3) = mapnode(i,j+1);
        elems(elem, 5) = 2;
        elem++;
    end
end
elem--; % number of elems



file = fopen('mixed_struct_grid.msh','w');
fprintf(file, '$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n%d\n',size(nodes,1));


for i = 1:size(nodes,1)
    fprintf(file, '%d %f %f %f\n', i, nodes(i,:));
end


fprintf(file, '$EndNodes\n$Elements\n');

len = size(lines,1) + size(elems,1);
fprintf(file, '%d\n', len);

% write lines and elems into "Elements"
for i = 1 : line
    fprintf(file, '%d %d %d %d %d %d %d\n', i, 1, 2, 0, lines(i,3), lines(i,1), lines(i,2));
    %              id  1  2  0 0/1 node1 node2
end

for i = 1 : elem
    if elems(i,5) == 3
        fprintf(file, '%d %d %d %d %d %d %d %d %d\n', i+line, 3, 2, 0, 0, elems(i,1), elems(i,2), elems(i,3), elems(i,4));
        %             id   3  2  0  0     nodes
    else
        fprintf(file, '%d %d %d %d %d %d %d %d\n', i+line, 2, 2, 0, 0, elems(i,1), elems(i,2), elems(i,3));
        %             id   3  2  0  0   nodes
    end
end

fprintf(file, '$EndElements');
fclose(file);


