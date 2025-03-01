// Parameters
r0 = 3;         // Radius to remove the singularity
r1 = 5;         // Radius of the intracellular circle
r2 = 6;        // Radius of the extracellular circle

// Characteristic length for mesh
lc = 0.1;

// Center of the circles
Point(1) = {0, 0, 0, lc};

// Points for the intracellular circle
Point(2) = {r1, 0, 0, lc};
Point(3) = {0, r1, 0, lc};
Point(4) = {-r1, 0, 0, lc};
Point(5) = {0, -r1, 0, lc};

// Points for the extracellular circle
Point(6) = {r2, 0, 0, lc};
Point(7) = {0, r2, 0, lc};
Point(8) = {-r2, 0, 0, lc};
Point(9) = {0, -r2, 0, lc};

// Points for the singularity circle
Point(10) = {r0, 0, 0, lc};
Point(11) = {0, r0, 0, lc};
Point(12) = {-r0, 0, 0, lc};
Point(13) = {0, -r0, 0, lc};

// Create the intracellular circle
Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};

// Create the extracellular circle
Circle(5) = {6, 1, 7};
Circle(6) = {7, 1, 8};
Circle(7) = {8, 1, 9};
Circle(8) = {9, 1, 6};

// Create the singularity circle
Circle(9) = {10, 1, 11};
Circle(10) = {11, 1, 12};
Circle(11) = {12, 1, 13};
Circle(12) = {13, 1, 10};

// Create the boundaries and surfaces
Curve Loop(1) = {1, 2, 3, 4};     // Boundary of the intracellular circle
Curve Loop(2) = {5, 6, 7, 8};     // Boundary of the extracellular circle
Curve Loop(3) = {9, 10, 11, 12};  // Boundary Singularity
Plane Surface(1) = {1, 3};          // Intracellular surface
Plane Surface(2) = {2, 1};       // Extracellular surface

// Label the boundaries
Physical Curve("cellMembrane", 3) = {1, 2, 3, 4};  // Boundary of the intracellular circle
Physical Curve("Boundary", 2) = {5, 6, 7, 8};      // Outer boundary of the extracellular circle
Physical Curve("InternalBoundaryCell", 4) = {9, 10, 11, 12}; // Internal boundary of the cell

// Label the surfaces
Physical Surface("intracelular", 3) = {1};
Physical Surface("extracelular", 1) = {2};
