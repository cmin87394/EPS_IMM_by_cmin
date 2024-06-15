
Jc = 0.04;   % Steering column moment of inertia
Kc = 15;   % Steering column stiffness
Bc = 0.072;   % Steering column viscous damping
N = 18.5;
Td = 0;   % Driver torque
Tf = 0;   % Friction torque

Jm = 0.0004;   % Motor moment of inertia
Bm = 0.0032;   % Motor shaft viscous damping
Kr = 81000;   % Tire spring rate
Mr = 32;   % Mass of the rack
Br = 3820;   % Viscous damping of the rack
Rp = 0.007;   % Steering column pinion radius
T = 0;   % Input torque of electrical power steering (EPS) systemT
Tr = 0;   % Road reaction torque on the rack and pinion
Teps = 0;   % Drive assistant torque

Jeq = Jm + (Rp^2 / N^2) * Mr;
Kn = (Kc + Kr * Rp^2) / N^2;
Beq = Bm + (Rp^2 / N^2) * Br;


a21 = -(Kc / Jc);
a22 = -(Bc / Jc);
a23 = Kc / (Jc * N);
d1 = (Td - Tf) / Jc;

a41 = Kc / (Jeq * N);
a43 = -(Kn / Jeq);
a44 = -(Beq / Jeq);
b4 = 1 / Jeq;
d2 = -(Rp / (Jeq * N)) * Tr + Teps;


A = [0 1 0 0; 
     a21 a22 a23 0; 
     0 0 0 1;
     a41 0 a43 a44];
B = [0 0 0 b4]';
C = [1 0 0 0];
D = 0;

sysc = ss(A,B,C,D);
sysd = c2d(sysc,0.01);
A2 = sysd.A;
B2 = sysd.B;
C2 = sysd.C;
D2 = sysd.D;

