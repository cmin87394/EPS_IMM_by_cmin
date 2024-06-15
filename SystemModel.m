function [A, B, C, D] = SystemModel()

run('states1')
run('states2')

A = {A1 A2};
B = {B1 B2};
C = {C1 C2};
D = {D1 D2};

end