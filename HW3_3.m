syms ux vx wx uy vy wy uz vz wz tx ty tz cp sp

A = [ux vx wx 0;
     uy vy wy 0;
     uz vz wz 0;
     0 0 0 1];

B = [ux vx wx;
     uy vy wy;
     uz vz wz];
 
C = [1 0 0 tx;
     0 1 0 ty;
     0 0 1 tz;
     0 0 0 1];

 disp(B \ B)
 
 
 A = [1 0 0;
      0 1 0;
      0 0 1];
 B = [0 0 ty;
      0 0 -tx;
      -ty tx 0];
 
%   disp(inv(B))

Rp = [cp -sp 0;
      sp cp 0;
      0 0 1];

  disp(Rp * B)