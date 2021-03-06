Point(1) = { -1,  1, 0 };
Point(2) = { -1, -1, 0 };
Point(3) = {  1, -1, 0 };
Point(4) = {  1,  0, 0 };
Point(5) = {  1,  1, 0 };
Point(6) = {  0,  1, 0 };
Line(1) = { 6, 1 };
Line(2) = { 1, 2 };
Line(3) = { 2, 3 };
Line(4) = { 3, 4 };
Circle(5) = { 4, 5, 6 };
Line Loop(7) = { 2, 3, 4, 5, 1 };
Plane Surface(7) = { 7 };
Physical Line("NeumannSW") = { 2, 3 };
Physical Line("Dirichlet") = { 1, 4 };
Physical Line("NeumannNE") = { 5 };
Physical Surface("Omega") = { 7 };
