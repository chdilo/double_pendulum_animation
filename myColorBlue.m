function h = myColorBlue(m)

n = fix(3/8*m);

b = [(1:n)'/n; ones(m-n,1)];
g = [zeros(n,1); (1:n)'/n; ones(m-2*n,1)];
r = [zeros(2*n,1); (1:m-2*n)'/(m-2*n)];

h = [r g b];

end