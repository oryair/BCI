
A = rand(100);
B = rand(100);

A = A * A';
B = B * B';

tic;
v1 = eig(B^-1 * A);
toc;

tic;
v2 = eig(A, B);
toc;

DDD = sort(v1) - sort(v2);
max(abs(DDD(:)))
