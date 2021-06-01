// Generate aquifer + cap mesh in gmsh with refinement
// vertically and horizontally around the screen interval
// Quad mesh within fine_r, tri mesh elsewhere.

//  5------line 3-----6-----line 12----9-----line 17----12    H_aq/2 + H_cap
//  |                 |                |                |
//  line 5          line 7            line 9         line 14
//  |                 |                |                |
//  32-----line 34----42               |                |     TOP CAP
//  |                 |                |                |
//  line 32         line 42            |                |
//  |                 |                |                |
//  3------line 2-----4-----line 11----8-----line 16----11    H_aq/2
//  |                 |                |                |
//  line 31         line 41            |                |
//  |                 |                |                |
//  31-----line 43----41               |                |     AQUIFER - TOP HALF
//  |                 |                |                |
//  line 4          line 6            line 8         line 13
//  |                 |                |                |
//  1------line 1-----2-----line 10----7-----line 15----10    zero
//  bh_r              fine_r           mid_r             max_r

// The above geometry is reflected in the plane y = 0 to generate the bottom half.

bh_r = 0.1;
mid_r = 50;
max_r = 1000;
//H_aq = 60;
H_cap = 40;

// Parameters to calculate thermal radius
cp_r = 800; // specific heat rock
cp_w = 4300; // specific heat water (representative value over temperature range)
//M_i = 1e6; // mass of injected fluid
rho_r = 2650; // density of rock
pi = 3.1416;
th_r = Sqrt(M_i*cp_w/(pi*H_aq*cp_r*rho_r)); // Thermal radius
fine_r = 2*th_r; // radius of fine mesh
Printf("Fine mesh radius is %f",fine_r);

// Mesh resolution in fine quad mesh areas
//ly_aq_medium = 0.3;
ly_aq_fine = 0.2;
ly_cap_fine = 0.2;
lr_fine = 0.3;
prog_y_top = 0.97;
prog_y_bottom = 1.0/prog_y_top;
//prog_r_top = 1.1;
prog_r_bottom = 1.0/prog_r_top;

// Height of fine area at top of aquifer and adjoining part of cap
fine_y_aq = 2;
fine_y_cap = 5;

ny_aq_medium = Ceil((H_aq/2 - fine_y_aq)/ly_aq_medium * 0.6);
ny_aq_fine = Ceil(fine_y_aq/ly_aq_fine);
ny_cap_fine = Ceil(fine_y_cap/ly_cap_fine * 0.8);
nr_fine = Ceil(fine_r/lr_fine);

top_aq = H_aq/2;
top = top_aq + H_cap;
top_aq_plus = top_aq + fine_y_cap;
top_aq_minus = top_aq - fine_y_aq;

// Characteristic lengths (in the quad mesh these are irrelevant)
l1 = 0.3;
l2 = 0.5;
l3 = 0.2;
l4 = 0.5;
l5 = 10;
l6 = 10;
l7 = 15;
l8 = 15;
l9 = 15;
l10 = 15;
l11 = 15;
l12 = 15;
l31 = l3;
l32 = 1;
l41 = l4;
l42 = 1;

Point(1) = {bh_r, 0, 0, l1};
Point(2) = {fine_r, 0, 0, l2};
Point(3) = {bh_r, top_aq, 0, l3};
Point(4) = {fine_r, top_aq, 0, l4};
Point(5) = {bh_r, top, 0, l5};
Point(6) = {fine_r, top, 0, l6};
Point(7) = {mid_r, 0, 0, l7};
Point(8) = {mid_r, top_aq, 0, l8};
Point(9) = {mid_r, top, 0, l9};
Point(10) = {max_r, 0, 0, l10};
Point(11) = {max_r, top_aq, 0, l11};
Point(12) = {max_r, top, 0, l12};
Point(31) = {bh_r, top_aq_minus, 0, l31};
Point(32) = {bh_r, top_aq_plus, 0, l32};
Point(41) = {fine_r, top_aq_minus, 0, l41};
Point(42) = {fine_r, top_aq_plus, 0, l42};


Line(1) = {1,2};
Line(2) = {3,4};
Line(3) = {5,6};
Line(4) = {1,31};
Line(5) = {32,5};
Line(6) = {2,41};
Line(7) = {42,6};
Line(8) = {7,8};
Line(9) = {8,9};
Line(10) = {2,7};
Line(11) = {4,8};
Line(12) = {6,9};
Line(13) = {10,11};
Line(14) = {11,12};
Line(15) = {7,10};
Line(16) = {8,11};
Line(17) = {9,12};
Line(31) = {31,3};
Line(32) = {3,32};
Line(34) = {32,42};
Line(41) = {41,4};
Line(42) = {4,42};
Line(43) = {31,41};

Line Loop(1) = {1,6,-43,-4}; // Aquifer fine lower
Line Loop(10) = {43,41,-2,-31}; // Aquifer fine upper
Line Loop(2) = {2,42,-34,-32}; // Cap fine lower
Line Loop(20) = {34,7,-3,-5}; // Cap fine upper
Line Loop(3) = {10,8,-11,-41,-6}; // Aquifer mid
Line Loop(4) = {11,9,-12,-7,-42}; // Cap mid
Line Loop(5) = {15,13,-16,-8}; // Aquifer coarse
Line Loop(6) = {16,14,-17,-9}; // Cap coarse

// Transfinite lines for quad mesh:
Transfinite Line{4,6} = ny_aq_medium Using Progression prog_y_top;
Transfinite Line{31,41} = ny_aq_fine;
Transfinite Line{32,42} = ny_cap_fine;
Transfinite Line{1,43,2,34} = nr_fine Using Progression prog_r_top;

// Surfaces with quad mesh:
Surface(101) = {1};
Surface(110) = {10};
Surface(102) = {2};
Transfinite Surface {101,110,102};
Recombine Surface {101,110,102};

// Surfaces with tri mesh:
Plane Surface(120) = {20};
Plane Surface(103) = {3};
Plane Surface(104) = {4};
Plane Surface(105) = {5};
Plane Surface(106) = {6};

Symmetry {0, 1, 0, 0} { Duplicata{ Surface{101,110,102,120,103,104,105,106}; } }

// Transfinite lines for quad mesh:
Transfinite Line{125,-123} = ny_aq_medium Using Progression prog_y_bottom;
Transfinite Line{130,128} = ny_aq_fine;
Transfinite Line{135,133} = ny_cap_fine;
Transfinite Line{124,129,134} = nr_fine Using Progression prog_r_bottom;

// Surfaces with quad mesh:
Transfinite Surface {121,126,131};
Recombine Surface {121,126,131};
