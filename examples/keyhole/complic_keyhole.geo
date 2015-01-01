Point(1) = {-1, 0, 0};
Point(2) = {-1,-2, 0};
Point(3) = { 1,-2, 0};
Point(4) = { 1, 0, 0};
Point(5) = { 0, 1, 0};
Point(6) = { 0, 1+Sqrt(2), 0};
Point(7) = { 0.5, 0.5, 0};
Point(8) = { 0, 1+Sqrt(2)/2, 0};
Point(9) = {-0.5, 0.5, 0};
Point(10) = { -0.5, -0.5, 0};
Point(11) = { -0.5, -1.5, 0};
Point(12) = { 0.5, -1.25, 0};
Point(13) = { 0.5, -0.75, 0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Circle(4) = {4, 5, 6};
Circle(5) = {6, 5, 1};
Circle(6) = {7, 5, 8};
Circle(7) = {8, 5, 9};
Circle(8) = {9, 5, 7};
Line(9) = { 10, 11 };
Line(10) = { 11, 12 };
Line(11) = { 12, 13 };
Line(12) = { 13, 10 };

Line Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
Plane Surface(1) = { 1 };
Line Loop(2) = { 9, 10, 11, 12 };
Plane Surface(2) = { 2 };

Physical Surface("LightBlue") = { 1 };
Physical Surface("LightGreen") = { 2 };
Physical Line("Violet") = { 1, 3 };
Physical Line("Black") = { 2 };
Physical Line("DarkBlue") = { 4, 5 };
Physical Line("Red") = { 6, 7, 8 };
