Mesh.MeshSizeFactor = 0.01025;

R1 = 1.13;
R2 = 1.83;
R3 = 2.05;

H1 = 0.94;
H2 = 1.13;

delta = 0.2;
R4 = R3 + delta;
H3 = H2 + delta;

Point(1) = {0, 0, 0, 1.0};
Point(2) = {0,  -H2, 0, 1.0};
Point(3) = {R3, -H2, 0, 1.0};
Point(4) = {R3,  H2, 0, 1.0};
Point(5) = {0,   H2, 0, 1.0};

Point(6) = {R1, -H1, 0, 1.0};
Point(7) = {R2, -H1, 0, 1.0};
Point(8) = {R1,  H1, 0, 1.0};
Point(9) = {R2,  H1, 0, 1.0};

Point(10) = {0,  -H3, 0, 1.0};
Point(11) = {R4, -H3, 0, 1.0};
Point(12) = {R4,  H3, 0, 1.0};
Point(13) = {0,   H3, 0, 1.0};

Point(14) = {0,   H1, 0, 1.0};
Point(15) = {0,   -H1, 0, 1.0};

Point(16) = {R1,   H2, 0, 1.0};
Point(17) = {R1,   -H2, 0, 1.0};

Point(18) = {R2,   H2, 0, 1.0};
Point(19) = {R2,   -H2, 0, 1.0};

Point(20) = {R3,   H1, 0, 1.0};
Point(21) = {R3,   -H1, 0, 1.0};

Point(22) = {R4,   H2, 0, 1.0};
Point(23) = {R4,   -H2, 0, 1.0};

Point(24) = {R3,   H3, 0, 1.0};
Point(25) = {R3,   -H3, 0, 1.0};

Point(26) = {R1,   H3, 0, 1.0};
Point(27) = {R1,   -H3, 0, 1.0};

Point(28) = {R2,   H3, 0, 1.0};
Point(29) = {R2,   -H3, 0, 1.0};

Point(30) = {R4,   H1, 0, 1.0};
Point(31) = {R4,   -H1, 0, 1.0};

Line(1) = {6, 7};
//+
Line(2) = {7, 9};
//+
Line(3) = {9, 8};
//+
Line(4) = {8, 6};
//+
Line(5) = {2, 17};
//+
Line(6) = {17, 19};
//+
Line(7) = {19, 3};
//+
Line(8) = {3, 21};
//+
Line(9) = {21, 20};
//+
Line(10) = {20, 4};
//+
Line(11) = {4, 18};
//+
Line(12) = {18, 16};
//+
Line(13) = {16, 5};
//+
Line(14) = {5, 14};
//+
Line(15) = {14, 15};
//+
Line(17) = {15, 2};
//+
Line(18) = {2, 10};
//+
Line(19) = {10, 27};
//+
Line(20) = {27, 29};
//+
Line(21) = {29, 25};
//+
Line(22) = {25, 11};
//+
Line(23) = {11, 23};
//+
Line(24) = {23, 31};
//+
Line(25) = {31, 30};
//+
Line(26) = {30, 22};
//+
Line(27) = {22, 12};
//+
Line(28) = {12, 24};
//+
Line(29) = {24, 28};
//+
Line(30) = {28, 26};
//+
Line(31) = {26, 13};
//+
Line(32) = {13, 5};
//+
Line(33) = {26, 16};
//+
Line(34) = {16, 8};
//+
Line(35) = {6, 17};
//+
Line(36) = {17, 27};
//+
Line(37) = {29, 19};
//+
Line(38) = {19, 7};
//+
Line(39) = {9, 18};
//+
Line(40) = {18, 28};
//+
Line(41) = {24, 4};
//+
Line(42) = {3, 25};
//+
Line(43) = {3, 23};
//+
Line(44) = {15, 6};
//+
Line(45) = {7, 21};
//+
Line(46) = {21, 31};
//+
Line(47) = {14, 8};
//+
Line(48) = {9, 20};
//+
Line(49) = {20, 30};
//+
Line(50) = {4, 22};


Curve Loop(1) = {47, 4, -44, -15};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {44, 35, -5, -17};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {1, -38, -6, -35};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {38, 45, -8, -7};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {45, 9, -48, -2};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {48, 10, 11, -39};
//+
Plane Surface(6) = {6};
//+
Curve Loop(7) = {12, 34, -3, 39};
//+
Plane Surface(7) = {7};
//+
Curve Loop(8) = {13, 14, 47, -34};
//+
Plane Surface(8) = {8};
//+
Curve Loop(9) = {5, 36, -19, -18};
//+
Plane Surface(9) = {9};
//+
Curve Loop(10) = {36, 20, 37, -6};
//+
Plane Surface(10) = {10};
//+
Curve Loop(11) = {21, -42, -7, -37};
//+
Plane Surface(11) = {11};
//+
Curve Loop(12) = {42, 22, 23, -43};
//+
Plane Surface(12) = {12};
//+
Curve Loop(13) = {43, 24, -46, -8};
//+
Plane Surface(13) = {13};
//+
Curve Loop(14) = {9, 49, -25, -46};
//+
Plane Surface(14) = {14};
//+
Curve Loop(15) = {49, 26, -50, -10};
//+
Plane Surface(15) = {15};
//+
Curve Loop(16) = {27, 28, 41, 50};
//+
Plane Surface(16) = {16};
//+
Curve Loop(17) = {29, -40, -11, -41};
//+
Plane Surface(17) = {17};
//+
Curve Loop(18) = {30, 33, -12, 40};
//+
Plane Surface(18) = {18};
//+
Curve Loop(19) = {31, 32, -13, -33};
//+
Plane Surface(19) = {19};



Transfinite Surface {1};
//+
Transfinite Surface {2};
//+
Transfinite Surface {9};
//+
Transfinite Surface {10};
//+
Transfinite Surface {3};
//+
Transfinite Surface {4};
//+
Transfinite Surface {11};
//+
Transfinite Surface {12};
//+
Transfinite Surface {13};
//+
Transfinite Surface {14};
//+
Transfinite Surface {5};
//+
Transfinite Surface {15};
//+
Transfinite Surface {16};
//+
Transfinite Surface {17};
//+
Transfinite Surface {6};
//+
Transfinite Surface {7};
//+
Transfinite Surface {18};
//+
Transfinite Surface {19};
//+
Transfinite Surface {8};

Physical Surface("core", 10) = {1, 8, 7, 6, 5, 4, 3, 2};
Physical Surface("reflector", 20) = {9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19};

Physical Curve("sym", 1) = {32, 14, 15, 17, 18};
Physical Curve("void", 2) = {31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19};
Physical Curve("blanket", 3) = {3, 4, 2, 1};
