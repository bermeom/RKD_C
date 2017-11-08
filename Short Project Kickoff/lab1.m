%% Lab1: ROBOTICS , KINEMATICS, DYNAMICS AND CONTROL
% Miguel Angle Bermeo Ayerbe
% install startup_rvc
% startup_rvc
close all;clear all;clc; 

%% Task 1 Drilling 40 holes
% Diameter of the drill 0.010 m. Number of holes: 20 in the front face and 20 in the rear face of half torus: 0.3m of tube diameter and 0.8 m outer radius . The pitch must be pi/19.
% 1.1 Load all torus points.

load('Drill_Mill_Shape.mat') % This is data.
figure(1);
image(Torus);
title('Torus');

%%
% Rotate all torus points using robotic toolbox


%% Task 2. Drilling and milling 25 shapes like the next figure
% 
scatter3(drill_mill_shape(1,:),drill_mill_shape(2,:),drill_mill_shape(3,:),'filled') % plotting the shape
axis equal

%% Help 1: Think about managing RTB data extructe.
% Drawing 30 points of circle at center point [2 3 4]', radius 5, using RTB data extructure
% 
n=30;px=2;py=3;pz=4;r=5;
for i=1:n
    Ptos(:,:,i)=transl(px,py,pz)*trotz((2*pi*i)/n)*transl(5,0,0);
end

coor_circle=transl(Ptos)';
figure();
scatter3(coor_circle(1,:),coor_circle(2,:),coor_circle(3,:),'r','LineWidth',2);

%% Help 3. What about a spiral in a cilinder
clear all; close all;clc; 
d=0.3;
n=8;
m=8*n;
r = 0.95;
r1 = 0.95-d/2;
r2 = d/2+0.001;
Ptos0 =[];
for i=0:(m)
    Ptos0(:,:,i+1)=troty((-pi*i)/m)*transl(r1,0,0)*trotz((2*pi*i)/n)*transl(r2,0,0);
end
spiral=transl(Ptos0)';
figure;

plot3(spiral(1,:),spiral(2,:),spiral(3,:),'g','LineWidth',2)
hold on
scatter3(spiral(1,:),spiral(2,:),spiral(3,:),'LineWidth',2);
% view (25,50)
xlabel('x');
ylabel('y');
zlabel('z');
grid on;
axis equal

% load('Data_groove_weld_fv_torus.mat') % This is some data you can use
% plot3(Weld_points(1,:),Weld_points(2,:),Weld_points(3,:),'r','LineWidth',2)
% hold on
% scatter3(Weld_points(1,:),Weld_points(2,:),Weld_points(3,:),'r','LineWidth',2)
% xlabel('x');
% ylabel('y');
% zlabel('z');
% axis 'equal'

m1 = 8;
Ptos1 = [];
for i=0:(m1-1)
    Ptos1(:,:,i*2+1)=troty((-pi*(i+1)+pi/2)/m1)*transl(r1,0,0)*trotz(-pi/2)*transl(r2,0,0);
    Ptos1(:,:,i*2+2)=troty((-pi*(i+1)+pi/2)/m1)*transl(r1,0,0)*trotz(pi/2)*transl(r2,0,0);
end
holes1=transl(Ptos1)';

% figure();
scatter3(holes1(1,:),holes1(2,:),holes1(3,:),'r','LineWidth',2);
xlabel('x');
ylabel('y');
zlabel('z');
grid on;
axis equal
% view (25,50)

m2 = 8;
n2 = 4;
Ptos=[];
for i=1:(m2)
        for (j=1:n2)
            Ptos(:,:,(i-1)*n2+j)=troty((-pi*(i))/(m2)+pi/24+(pi*j)/100)*transl(r,0,0);%*trotz(2*pi)*transl(r2,0,0);
        end
end
holes2=transl(Ptos)';

% figure();
scatter3(holes2(1,:),holes2(2,:),holes2(3,:),'m','LineWidth',2);
xlabel('x');
ylabel('y');
zlabel('z');
grid on;
axis equal

fv=stlread('Torus.stl');% fv is a struct with faces and vertices
fh = [fv.vertices'; ones(1,size(fv.vertices,1))];
fh=troty(pi)*fh;
fh(1:3,:) = (r).*(fh(1:3,:)./(max(max(fh(1:3,:)))));
fh = transl(r1,0,0)*fh;
fv.vertices = fh(1:3,:)';
% figure();
patch(fv,'FaceColor',       [0.7 0.7 1.0], ...
         'EdgeColor',       'none',        ...
         'FaceLighting',    'gouraud',     ...
         'AmbientStrength', 0.08);
% Add a camera light, and tone down the specular highlighting
camlight('headlight');
material('dull');
hold on;
grid on;
% Fix the axes scaling, and set a nice view angle
axis('image');
axis 'equal';
% alpha 0.7;
view([-135 35]);


%%
figure
xlabel('x');
ylabel('y');
zlabel('z');
axis 'equal';
table = [1 3 3 1;0 0 5 5;4 4 6 6] 
fill3(table(1,:),table(2,:),table(3,:),'r');
alpha 0.3;

norm(sum(table(1,:))/3)
theta = pi*(30/180);%% 30
fv=stlread('Torus.stl');% fv is a struct with faces and vertices
fh = [fv.vertices'; ones(1,size(fv.vertices,1))];
fh=troty(pi)*fh;
fh(1:3,:) = (r).*(fh(1:3,:)./(max(max(fh(1:3,:)))));
fh = transl(r/2+sum(table(1,:)-min(table(1,:)))/3+min(table(1,:)),sum(table(2,:)-min(table(2,:)))/3+min(table(2,:)),sum(table(3,:)-min(table(3,:)))/3+min(table(3,:)))*trotx(theta)*fh;
fv.vertices = fh(1:3,:)';
% figure();
patch(fv,'FaceColor',       [0.7 0.7 1.0], ...
         'EdgeColor',       'none',        ...
         'FaceLighting',    'gouraud',     ...
         'AmbientStrength', 0.08);
% Add a camera light, and tone down the specular highlighting
camlight('headlight');
material('dull');
hold on;
grid on;
% Fix the axes scaling, and set a nice view angle
axis('image');
axis 'equal';
% alpha 0.7;
view([-135 35]);





%% Animating the Puma manipulator performing a task
% *Initialize the workspace*:
%%
% We clear the workspace and set the value of some constants of the problem
% clc;clear;close all;
radius=0.20;
n=200;
INI = transl(-0.25, 0.25,-0.5); %center of the part
% Plot a sketch of the environment of the robot
mdl_puma560
% figure(6)
p560.plot(qz);
hold on;
circle1 = holes1;%%circle([-0.25 0.25 -0.5], radius);
circle2 = circle([-0.25 0.25 -0.5], 0.5*radius);
plot3(circle1(1,:), circle1(2,:), circle1(3,:),'g','LineWidth',1);
patch(circle2(1,:), circle2(2,:), circle2(3,:),'r');

for i=1:n
    Laser_Pose(:,:,i)= INI*trotx(-pi/2)*troty(2*pi*i/n)*transl(0, 0, -radius);
end

Q= p560.ikine6s(INI*holes1, 'run');
p560.plot(Q);


