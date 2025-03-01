// Create cuboid of idealized cells

// Characteristic length for mesh
lca = 0.14;

// Adjust to change the origin of the mesh
x0=-1;
y0=-0.5;
z0=-0.5;

Nx = 2;
Ny = 2;
Nz = 1;

// Dimensions of cell
h = 1;
w = 1;
l = 1;

// Connection in x direction
hx = 1/2;
wx = 1/2;
lx = 1/2;

// Connection in y direction
hy = 1/2;
wy= 1/2;
ly = 1/2;

// Connection in z direction
hz= 1/2;
wz = 1/2;
lz = 1/2;

i2 = 0;
Macro Cell
 i = 0;
 If (nx == 1)
 pa = newp; Point(pa) = {x+x0, y+y0+w/2+wx/2, z+z0+h/2-hx/2, lca};
 pb = newp; Point(pb) = {x+x0, y+y0+w/2-wx/2, z+z0+h/2-hx/2, lca};
 pc = newp; Point(pc) = {x+x0, y+y0+w/2-wx/2, z+z0+h/2+hx/2, lca};
 pd = newp; Point(pd) = {x+x0, y+y0+w/2+wx/2, z+z0+h/2+hx/2, lca};
 l1 = newl; Line(l1) = {pa, pb};
 l2 = newl; Line(l2) = {pb, pc};
 l3 = newl; Line(l3) = {pc, pd};
 l4 = newl; Line(l4) = {pd, pa};
 c = newc; Curve Loop(c) = {l1, l2, l3, l4};
 s = news; Plane Surface(s) = {c};
 membrane~{t}[i] = s;
 membrane[i2] = s;
 i += 1;
 i2 += 1;
 EndIf
 //If (nx != 1)
    out[] = Extrude{lx, 0, 0}{Curve{l1}; Curve{l2};Curve{l3}; Curve{l4};};
 For j In {1:4}
  membrane~{t}[i] = out[-3+4*j];
  membrane[i2] = out[-3 + 4*j];
  i2 += 1;
  i += 1;
 EndFor
 //EndIf
 p1 = newp; Point(p1) = {x+x0+lx,y+y0,z+z0,lca};
 p2 = newp; Point(p2) = {x+x0+lx, y+y0+w, z+z0, lca};
 p3 = newp; Point(p3) = {x+x0+lx, y+y0+w, z+z0+h, lca};
 p4 = newp; Point(p4) = {x+x0+lx, y+y0, z+z0+h, lca};
 l1 = newl; Line(l1) = {p1, p2};
 l2 = newl; Line(l2) = {p2, p3};
 l3 = newl; Line(l3) = {p3, p4};
 l4 = newl; Line(l4) = {p4, p1};
 c2 = newc; Curve Loop(c2) = {l1, l2, l3, l4};
 //If (nx != 1)
 c3 = newc; Curve Loop(c3) = {out[0], out[4], out[8], out[12]};
 s = news; Plane Surface(s) = {c2, c3};
 //Else
 //s = news; Plane Surface(s) = {c2};
 //EndIf
 membrane~{t}[i] = s;
 membrane[i2] = s;
 i2 += 1;
 i += 1;
 out[] = Extrude{l, 0, 0}{Curve{l1}; Curve{l2};Curve{l3}; Curve{l4};};
 For j In {1:4}
  membrane~{t}[i] = out[-3+4*j];
  i += 1;
  membrane[i2] = out[-3+4*j];
  i2 += 1;
 EndFor
 
 //If (ny != 1)
 Delete{ Surface{out[13]};}
 If (ny != 1)
 outy2[] = Extrude{0, wy, 0}{Curve{outy~{nx}[0]}; Curve{outy~{nx}[4]}; Curve{outy~{nx}[8]}; Curve{outy~{nx}[12]};};
 Else
 pya = newp; Point(pya) = {x+x0 + lx + l/2 - ly/2, y+y0-wy, z+z0 + h/2 - hy/2, lca};
 pyb = newp; Point(pyb) = {x+x0 + lx + l/2 - ly/2, y+y0-wy, z+z0+h/2+hy/2, lca};
 pyc = newp; Point(pyc) = {x+x0 + lx + l/2 + ly/2, y+y0-wy, z+z0+h/2+hy/2, lca};
 pyd = newp; Point(pyd) = {x+x0 + lx + l/2 + ly/2, y+y0-wy, z+z0+h/2-hy/2, lca};
 l1y = newl; Line(l1y) = {pya, pyb};
 l2y = newl; Line(l2y) = {pyb, pyc};
 l3y = newl; Line(l3y) = {pyc, pyd};
 l4y = newl; Line(l4y) = {pyd, pya};
 outy2[] = Extrude{0, wy, 0}{Curve{l1y}; Curve{l2y}; Curve{l3y}; Curve{l4y};};
 EndIf
 For j In {1:4}
 membrane~{t}[i] = outy2[-3+4*j];
 membrane[i2] = outy2[-3+4*j];
 i += 1;
 i2 += 1;
 EndFor
 c2 = newc; Curve Loop(c2) = {l4, out[14], -out[12], out[15]};
 c3 = newc; Curve Loop(c3) = {outy2[0], outy2[4], outy2[8], outy2[12]};
 Plane Surface(out[13]) = {c2, c3};
 If (ny != 1)
 gapy_s = gapy[nx];
 Else
 cy = newc; Curve Loop(cy) = {l1y, l2y, l3y, l4y};
 gapy_s = news; Plane Surface(gapy_s) = {cy};
 EndIf
 //If (ny != Ny)
 Delete{ Surface{out[5]};}
 pa = newp; Point(pa) = {x+x0 + lx + l/2 - ly/2, y+y0+w, z+z0+h/2-hy/2,lca};
 pb = newp; Point(pb) = {x+x0 + lx + l/2 - ly/2, y+y0+w, z+z0+h/2+hy/2, lca};
 pc = newp; Point(pc) = {x+x0 + lx + l/2 + ly/2, y+y0+w, z+z0+h/2+hy/2, lca};
 pd = newp; Point(pd) = {x+x0 + lx + l/2 + ly/2, y+y0+w, z+z0+h/2-hy/2, lca};
 la = newl; Line(la) = {pa, pb};
 lb = newl; Line(lb) = {pb, pc};
 lc = newl; Line(lc) = {pc, pd};
 ld = newl; Line(ld) = {pd, pa};
 c2 = newc; Curve Loop(c2) = {l2, out[6], -out[4], out[7]};
 c3 = newc; Curve Loop(c3) = {la, lb, lc, ld};
 Plane Surface(out[5]) = {c2, c3};
 outy[] = Extrude{0,wy,0}{Curve{la}; Curve{lb}; Curve{lc}; Curve{ld};};
 c = newc; Curve Loop(c) = {outy[0], outy[4], outy[8], outy[12]};
 outy~{nx}[] = outy[];
 For j In {1:4}
 membrane~{t}[i] = outy[-3+4*j];
 membrane[i2] = outy[-3+4*j];
 i += 1;
 i2 += 1;
 EndFor
 gapy[nx] = news; Plane Surface(gapy[nx]) = {c};
 If (ny != Ny)
 Physical Surface(gap_i) = {gapy[nx]};
 gap_i += 1;
 Else
 membrane~{t}[i] = gapy[nx];
 membrane[i2] = gapy[nx];
 i2 += 1;
 i += 1;
 EndIf
 //EndIf
 If (nz != 1)
 Delete{ Surface{out[1]};}
 outz2[] = Extrude{0,0,hz}{Curve{outz~{ny}~{nx}[0]};Curve{outz~{ny}~{nx}[4]};Curve{outz~{ny}~{nx}[8]};Curve{outz~{ny}~{nx}[12]};};
 For j In {1:4}
 membrane~{t}[i] = outz2[-3+4*j];
 membrane[i2] = outz2[-3+4*j];
 i += 1;
 i2 += 1;
 EndFor
 c2 = newc; Curve Loop(c2) = {l1, out[2], -out[0], out[3]};
 c3 = newc; Curve Loop(c3) = {outz2[0], outz2[4], outz2[8], outz2[12]};
 Plane Surface(out[1]) = {c2, c3};
 gapz_s = gapz~{ny}[nx];
 EndIf
 If (nz != Nz)
 Delete{ Surface{out[9]};}
 pa = newp; Point(pa) = {x+x0 + lx + l/2 - lz/2, y+y0 + w/2 -wz/2, z+z0+h, lca};
 pb = newp; Point(pb) = {x+x0 + lx + l/2 - lz/2, y+y0 + w/2 + wz/2, z+z0+h, lca};
 pc = newp; Point(pc) = {x+x0 + lx + l/2 + lz/2, y+y0 + w/2 + wz/2, z+z0+h, lca};
 pd = newp; Point(pd) = {x+x0 + lx + l/2 + lz/2, y+y0 + w/2 - wz/2, z+z0 + h, lca};
 la = newl; Line(la) = {pa, pb};
 lb = newl; Line(lb) = {pb, pc};
 lc = newl; Line(lc) = {pc, pd};
 ld = newl; Line(ld) = {pd, pa};
 c2 = newc; Curve Loop(c2) = {l3, out[10], -out[8], out[11]};
 c3 = newc; Curve Loop(c3) = {la, lb, lc, ld};
 Plane Surface(out[9]) = {c2, c3};
 outz[] = Extrude{0,0,hz}{Curve{la}; Curve{lb}; Curve{lc}; Curve{ld};};
 c = newc; Curve Loop(c) = {outz[0], outz[4], outz[8], outz[12]};
 outz~{ny}~{nx}[] = outz[];
 For j In {1:4}
 membrane~{t}[i] = outz[-3+4*j];
 membrane[i2] = outz[-3+4*j];
 i += 1;
 i2 += 1;
 EndFor
 gapz~{ny}[nx] = news; Plane Surface(gapz~{ny}[nx]) = {c};
 Physical Surface(gap_i) = {gapz~{ny}[nx]};
 gap_i += 1;
 EndIf
 //If (nx != Nx)

 
 pa = newp; Point(pa) = {x+x0+l+lx, y+y0+w/2+wx/2, z+z0+h/2-hx/2, lca};
 pb = newp; Point(pb) = {x+x0+l+lx, y+y0+w/2-wx/2, z+z0+h/2-hx/2, lca};
 pc = newp; Point(pc) = {x+x0+l+lx, y+y0+w/2-wx/2, z+z0+h/2+hx/2, lca};
 pd = newp; Point(pd) = {x+x0+l+lx, y+y0+w/2+wx/2, z+z0+h/2+hx/2, lca};
 la = newl; Line(la) = {pa, pb};
 lb = newl; Line(lb) = {pb, pc};
 lc = newl; Line(lc) = {pc, pd};
 ld = newl; Line(ld) = {pd, pa};
 out1[] = Extrude{lx, 0, 0} {Curve{la}; Curve{lb}; Curve{lc}; Curve{ld};};
 For j In {1:4}
  membrane~{t}[i] = out1[-3+4*j];
  i += 1;
  membrane[i2] = out1[-3+4*j];
  i2 += 1;
 EndFor
 ca = newc; Curve Loop(ca) = {out[0], out[4], out[8], out[12]};
 cb = newc; Curve Loop(cb) = {la, lb, lc, ld};
 s = news; Plane Surface(s) = {ca, cb};
 //Else
 //ca = newc; Curve Loop(ca) = {out[0], out[4], out[8], out[12]};
 //s = news; Plane Surface(s) = {ca};
 //EndIf
 membrane~{t}[i] = s;
 membrane[i2] = s;
 i2 += 1;
 i += 1;
 If (nx != Nx)
 l1 = out1[0];
 l2 = out1[4];
 l3 = out1[8];
 l4 = out1[12];
 g = newc; Curve Loop(g) = {out1[0], out1[4], out1[8], out1[12]};
 gap2 = news; Plane Surface(gap2) = {g};
   Physical Surface(gap_i) = {gap2};
   gap_i += 1;
 Else
 l1 = out1[0];
 l2 = out1[4];
 l3 = out1[8];
 l4 = out1[12];
 g = newc; Curve Loop(g) = {out1[0], out1[4], out1[8], out1[12]};
 gap2 = news; Plane Surface(gap2) = {g};
 membrane~{t}[i] = gap2;
 membrane[i2] = gap2;
 i2 += 1;
 i += 1;
 EndIf
 If (ny == 1)
 membrane~{t}[i] = gapy_s;
 membrane[i2] = gapy_s;
 i2 += 1;
 i += 1;
 EndIf
 Physical Surface(t+2) = membrane~{t}[];
 If (ny != Ny)
   membrane~{t}[i] = gapy[nx];
   i += 1;
 EndIf
 If (ny != 1)
 membrane~{t}[i] = gapy_s;
 i += 1;
   EndIf
 If (nz != Nz)
 membrane~{t}[i] = gapz~{ny}[nx];
 i+=1;
 EndIf
 If (nz != 1)
 membrane~{t}[i] = gapz_s;
 i+=1;
 EndIf
 If (nx != Nx)
   membrane~{t}[i] = gap2;
   i += 1;
 EndIf
 If (nx != 1)
   membrane~{t}[i] = gap;
   i += 1;
 EndIf
 s = news; Surface Loop(t+2) = membrane~{t}[];
 Volume(t+2) = {t+2};
 Physical Volume(t+2) = {t+2};
 
 gap = gap2;
Return

x = 0; y=0; z = 0;

//Nx = 3;
gap_i = Nx*Ny*Nz + 3;
//pa = newp; Point(pa) = {x-wx, y + h/2 - hx/2, 0, lca};
//pb = newp; Point(pb) = {x-wx, y + h/2 + hx/2, 0, lca};
//base = newl; Line(base) = {pa, pb};
For nz In {1:Nz}
For ny In {1:Ny}
For nx In {1:Nx}
  t = (nz-1) * (Nx * Ny) + (ny - 1) * Nx + nx;
  Call Cell;
  x += l + lx + lx;

EndFor
y += w + wy + wy;
x = 0;
EndFor
y=0;
x = 0;
z += h + hz + hz;
EndFor

ec_x = 3/8;
ec_y = 3/8;
ec_z = 3/8;

p1 = newp; Point(p1) = {-ec_x+x0, -ec_y-wy+y0, -ec_z+z0, lca};
p2 = newp; Point(p2) = {-ec_x+x0, (w + 2 * wy) * Ny + ec_y - wy+y0, -ec_z+z0, lca};
p3 = newp; Point(p3) = {- ec_x+x0, (w + 2 * wy) * Ny +ec_y -wy+y0, (h + 2 * hz) * Nz - 2*hz+ec_z+z0, lca};
p4 = newp; Point(p4) = {-ec_x+x0, -ec_y-wy+y0, (h + 2*hz) * Nz - 2*hz+ec_z+z0, lca};
l1 = newl; Line(l1) = {p1, p2};
l2 = newl; Line(l2) = {p2, p3};
l3 = newl; Line(l3) = {p3, p4};
l4 = newl; Line(l4) = {p4, p1};
out[]=Extrude{Nx * (l + 2 * lx) + 2* ec_x, 0, 0} {Curve{l1}; Curve{l2};Curve{l3}; Curve{l4};};

c = newc; Curve Loop(c) = {l1, l2, l3, l4};
s = news; Plane Surface(s) = {c};

c = newc; Curve Loop(c) = {out[0], out[4], out[8], out[12]};
s2 = news; Plane Surface(s2) = {c};

s3 = news; Surface Loop(1) = {s, s2, out[1], out[5], out[9], out[13]};
Surface Loop(s3) = membrane[];

// Here we define boundary regions for Dirichlet boundary conditions (region 1) and Neummann (region 2)
Physical Surface(2) = {out[1], out[5], out[9], out[13]};
Physical Surface(1) = {s, s2};
Volume(1) = {1, s3};
Physical Volume(1) = {1};

Point(0123456789) = {0, 0, 0, lca};
