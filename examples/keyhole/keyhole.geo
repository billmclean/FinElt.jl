Point(1) = {-1, 0, 0};
Point(2) = {-1,-2, 0};
Point(3) = { 1,-2, 0};
Point(4) = { 1, 0, 0};
Point(5) = { 0, 1, 0};
Point(6) = { 0, 1+Sqrt(2), 0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Circle(4) = {4, 5, 6};
Circle(5) = {6, 5, 1};

Line Loop(7) = {1, 2, 3, 4, 5};
Plane Surface(1) = { 7 };

Physical Surface("Omega") = { 1 };
Physical Line("Gamma") = { 1, 2, 3, 4, 5 };
