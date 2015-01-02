Lx = 2.0;
Ly = 1.0;
Point(1) = { 0.0, 0.0, 0.0 };
Point(2) = { Lx,  0.0, 0.0 };
Point(3) = { Lx,  Ly,  0.0 };
Point(4) = { 0.0, Ly,  0.0 };

Line(1) = { 1, 2 };
Line(2) = { 2, 3 };
Line(3) = { 3, 4 };
Line(4) = { 4, 1 };

Line Loop(1) = { 1, 2, 3, 4 };
Plane Surface(1) = { 1 };

Physical Surface("Omega") = { 1 };

Physical Line("Top")    = 3;
Physical Line("Bottom") = 1;
Physical Line("Left")   = 4;
Physical Line("Right")  = 2;
