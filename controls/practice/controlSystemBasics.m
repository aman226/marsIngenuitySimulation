% m                   mass                                    1.0 kg
% k                   spring constant                         1.0 N/m
% b                   damping constant                        0.2 Ns/m
% F                   input force                             1.0 N
%%
m = 1.00;
k = 1.00;
b = 0.200;
F = 1.00;

A = [0 1.0; -k/m -b/m];
B = [0 1/m]';
C = [1.0 0];
D = 0;
sys = ss(A,B,C,D);

%%
rlocus(sys);

num = 1.0;
den = [m b k];
sys = tf(num,den);

s = tf('s');
ans = ((s*eye(size(A,1)) - A)*[5/s;6/s])

pole(sys);
eig(A);
%%