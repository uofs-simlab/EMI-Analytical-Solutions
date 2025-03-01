// Parameters
r0 = 3;         // Radius to remove the singularity
r1 = 5;         // Radius of the intracellular sphere
r2 = 6;        // Radius of the extracellular sphere

// Characteristic length for mesh
lc = 0.1;

// Center of the spheres
Point(1) = {0, 0, 0, lc};

// Points for the middle sphere shell (radius r1)
Point(2) = {r1, 0, 0, lc};
Point(3) = {0, r1, 0, lc};
Point(4) = {-r1, 0, 0, lc};
Point(5) = {0, -r1, 0, lc};
Point(6) = {0, 0, r1, lc};
Point(7) = {0, 0, -r1, lc};

// Points for the inner sphere shell (radius r0)
Point(8) = {r0, 0, 0, lc};
Point(9) = {0, r0, 0, lc};
Point(10) = {-r0, 0, 0, lc};
Point(11) = {0, -r0, 0, lc};
Point(12) = {0, 0, r0, lc};
Point(13) = {0, 0, -r0, lc};

// Points for the outer sphere shell (radius r2)
Point(14) = {r2, 0, 0, lc};
Point(15) = {0, r2, 0, lc};
Point(16) = {-r2, 0, 0, lc};
Point(17) = {0, -r2, 0, lc};
Point(18) = {0, 0, r2, lc};
Point(19) = {0, 0, -r2, lc};

// Circles for the middle sphere shell
Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};
Circle(5) = {5, 1, 7};
Circle(6) = {7, 1, 3};
Circle(7) = {3, 1, 6};
Circle(8) = {6, 1, 5};
Circle(9) = {4, 1, 7};
Circle(10) = {7, 1, 2};
Circle(11) = {2, 1, 6};
Circle(12) = {6, 1, 4};

// Circles for the inner sphere shell
Circle(13) = {8, 1, 9};
Circle(14) = {9, 1, 10};
Circle(15) = {10, 1, 11};
Circle(16) = {11, 1, 8};
Circle(17) = {11, 1, 13};
Circle(18) = {13, 1, 9};
Circle(19) = {9, 1, 12};
Circle(20) = {12, 1, 11};
Circle(21) = {10, 1, 13};
Circle(22) = {13, 1, 8};
Circle(23) = {8, 1, 12};
Circle(24) = {12, 1, 10};

// Circles for the outer sphere shell
Circle(25) = {14, 1, 15};
Circle(26) = {15, 1, 16};
Circle(27) = {16, 1, 17};
Circle(28) = {17, 1, 14};
Circle(29) = {17, 1, 19};
Circle(30) = {19, 1, 15};
Circle(31) = {15, 1, 18};
Circle(32) = {18, 1, 17};
Circle(33) = {16, 1, 19};
Circle(34) = {19, 1, 14};
Circle(35) = {14, 1, 18};
Circle(36) = {18, 1, 16};

// Create surfaces for the middle sphere shell
Line Loop(1) = {2, 9, 6};
Ruled Surface(1) = {1};
Line Loop(2) = {12, -2, 7};
Ruled Surface(2) = {2};
Line Loop(3) = {1, -6, 10};
Ruled Surface(3) = {3};
Line Loop(4) = {10, -4, 5};
Ruled Surface(4) = {4};
Line Loop(5) = {1, 7, -11};
Ruled Surface(5) = {5};
Line Loop(6) = {11, 8, 4};
Ruled Surface(6) = {6};
Line Loop(7) = {5, -9, 3};
Ruled Surface(7) = {7};
Line Loop(8) = {8, -3, -12};
Ruled Surface(8) = {8};

// Create surfaces for the inner sphere shell
Line Loop(9) = {14, 21, 18};
Ruled Surface(9) = {9};
Line Loop(10) = {24, -14, 19};
Ruled Surface(10) = {10};
Line Loop(11) = {13, -18, 22};
Ruled Surface(11) = {11};
Line Loop(12) = {22, -16, 17};
Ruled Surface(12) = {12};
Line Loop(13) = {13, 19, -23};
Ruled Surface(13) = {13};
Line Loop(14) = {23, 20, 16};
Ruled Surface(14) = {14};
Line Loop(15) = {17, -21, 15};
Ruled Surface(15) = {15};
Line Loop(16) = {20, -15, -24};
Ruled Surface(16) = {16};

// Create surfaces for the outer sphere shell
Line Loop(17) = {26, 33, 30};
Ruled Surface(17) = {17};
Line Loop(18) = {36, -26, 31};
Ruled Surface(18) = {18};
Line Loop(19) = {25, -30, 34};
Ruled Surface(19) = {19};
Line Loop(20) = {34, -28, 29};
Ruled Surface(20) = {20};
Line Loop(21) = {25, 31, -35};
Ruled Surface(21) = {21};
Line Loop(22) = {35, 32, 28};
Ruled Surface(22) = {22};
Line Loop(23) = {29, -33, 27};
Ruled Surface(23) = {23};
Line Loop(24) = {32, -27, -36};
Ruled Surface(24) = {24};

// Gap junction 
Curve Loop(25) = {8, 5, 6, 7};
Curve Loop(26) = {20, 17, 18, 19};
Plane Surface(27) = {25, 26};
Physical Curve("Gap", 8) = {20, 17, 18, 19};
Physical Surface("Gap", 5) = {27}; 


// Cells
Surface Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8}; //Both cells
Surface Loop(3) = {1, 2, 8, 7}; // Cell1
Surface Loop(4) = {5, 6, 4, 3}; // Cell2
Physical Surface("cellMembrane1", 3) = {1, 2, 7, 8};
Physical Surface("cellMembrane2", 4) = {5, 6, 4, 3};

Surface Loop(5) = {9, 10, 11, 12, 13, 14, 15, 16}; // Both InternalBoundaryCells
Surface Loop(6) = {9, 10, 15, 16}; // InternalBoundaryCell1
Surface Loop(7) = {13, 14, 11, 12}; // InternalBoundaryCell2
Physical Surface("internalBoundaryCell1", 6) = {9, 10, 15, 16};
Physical Surface("internalBoundaryCell2", 7) = {13, 14, 11, 12};

Surface Loop(8) = {1, 2, 8, 7, 9, 10, 15, 16, 27}; 
Volume(1) = {8}; // intracelullarCell1
Physical Volume("intracelullarCell1", 3) = {1};

Surface Loop(9) = {5, 6, 4, 3, 13, 14, 11, 12, 27}; 
Volume(2) = {9}; // intracelullarCell2
Physical Volume("intracelullarCell2", 4) = {2};


//Extracelullar space
Surface Loop(10) = {1, 2, 3, 4, 5, 6, 7, 8}; // Cell membranes
Surface Loop(11) = {17, 18, 19, 20, 21, 22, 23, 24}; // Boundary
Volume(3) = {11, 10}; // Outer volume between outer shell and middle shell
Physical Volume("extracellular", 1) = {3};
Physical Surface("boundary", 2) = {17, 18, 19, 20, 21, 22, 23, 24};

