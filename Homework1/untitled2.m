%% Homework #1
% Miguel Angel Bermeo Ayerbe
%%
% Given a transformation matrix:
A_B_T = [  2^(-1/2) 0        2^(-1/2) 1 ;...
          -1/2     2^(-1/2)  1/2      2 ;...
          -1/2    -2^(-1/2)  1/2      3 ;...
           0       0         0        1  ]

%% 1. Shows that $_{A}^{B}\textrm{R}_{3x3}$ is a rotation matrix
% Since its the three columns are orthogonal, then a consequence of this is that:
%%
% $_{A}^{B}\textrm{R}_{}^{-1} = _{A}^{B}\textrm{R}_{}^{T}$  and  $det(_{A}^{B}\textrm{R}) = 1$ 

A_B_R = t2r(A_B_T)
A_B_R*A_B_R'
det(A_B_R)


%% 2. What is the meaning of columns and rows of  $_{A}^{B}\textrm{R}_{3x3}$ and $_{}^{A}\textrm{M}$
% $_{A}^{B}\textrm{R}_{3x3}$ :
% Each column is a vector in the order X, Y and Z, and the components of any vector are simply the projections of that vector onto the unit directions of its reference frame
% $_{}^{A}\textrm{M}$:
% This means that the components of ${P}^{A}$ have numerical values that indicate distances along the
% axes of {A}

%% 3.

B_P = [ 4 5 6 1 ]';
A_P = A_B_T*B_P

%% 4.

A_P_1 = [ 4 5 6 1 ]';
T = A_B_T;
A_P_2 = T*A_P_1

%% 5. 
B_A_T = inv(A_B_T)

%% 6. 

betha = atan2(-A_B_R(3,1),sqrt(A_B_R(1,1)^2+A_B_R(2,1)^2))
alpha = atan2( A_B_R(2,1)/cos(betha),A_B_R(1,1)/cos(betha))
gamma = atan2( A_B_R(3,2)/cos(betha),A_B_R(3,3)/cos(betha))
my_xyz_fangle = [gamma betha alpha]
xyz_fangle= tr2rpy(A_B_R)
R = rpy2r(xyz_fangle(1), xyz_fangle(2), xyz_fangle(3))
A_B_R

[theta, v] = tr2angvec(A_B_R)


%% 7. 

A_B_r1 = A_B_R(:,1)
A_B_r2 = A_B_R(:,2)
A_B_r3 = A_B_R(:,3)

norm(A_B_r1)
norm(A_B_r2)
norm(A_B_r3)

%[vectors,values] = eig(A_B_R)
[theta, v] = tr2angvec(A_B_R)

%% 8.
gamma = tr2eul(A_B_R)