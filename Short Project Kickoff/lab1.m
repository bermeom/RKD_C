%% Lab1: ROBOTICS , KINEMATICS, DYNAMICS AND CONTROL
% Miguel Angle Bermeo Ayerbe
% install startup_rvc
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
n=8;px=2;py=3;pz=4;r1=0.95;r2=d/2;h=9;
m=100;
k = 7;
for i=0:(m)
    Ptos(:,:,i+1)=troty((-pi*i)/m)*transl(r1,0,0)*trotz((2*pi*i)/n-pi)*transl(r2,0,0);
end
coor_circle=transl(Ptos)';
figure;

plot3(coor_circle(1,:),coor_circle(2,:),coor_circle(3,:),'g','LineWidth',2)
hold on
scatter3(coor_circle(1,:),coor_circle(2,:),coor_circle(3,:),'LineWidth',2);
% view (25,50)
xlabel('x');
ylabel('y');
zlabel('z');
grid on;
axis equal

%% Reading '.stl' file
fv=stlread('Torus.stl')% fv is a struct with faces and vertices
figure(4)
patch(fv,'FaceColor',       [0.8 0.8 1.0], ...
         'EdgeColor',       'none',        ...
         'FaceLighting',    'gouraud',     ...
         'AmbientStrength', 0.15);
% Add a camera light, and tone down the specular highlighting
camlight('headlight');
material('dull');
% Fix the axes scaling, and set a nice view angle
axis('image');
axis 'equal'
view([-135 35]);

fh = [fv.vertices'; ones(1,size(fv.vertices,1))];
fh1=troty(pi)*fh;
fv.vertices = fh1(1:3,:)';
figure(5);
patch(fv,'FaceColor',       [0.8 0.8 1.0], ...
         'EdgeColor',       'none',        ...
         'FaceLighting',    'gouraud',     ...
         'AmbientStrength', 0.15);
% Add a camera light, and tone down the specular highlighting
camlight('headlight');
material('dull');

% Fix the axes scaling, and set a nice view angle
axis('image');
axis 'equal'
view([-135 35]);

%% Animating the Puma manipulator performing a task
% *Initialize the workspace*:
%%
% We clear the workspace and set the value of some constants of the problem
clc;
clear;
close all;
radius=0.20;
n=200;
INI = transl(-0.25, 0.25,-0.5); %center of the part
% Plot a sketch of the environment of the robot
mdl_puma560
figure(6)
p560.plot(qz);
hold on;
circle1 = circle([-0.25 0.25 -0.5], radius);
circle2 = circle([-0.25 0.25 -0.5], 0.5*radius);
plot3(circle1(1,:), circle1(2,:), circle1(3,:),'g','LineWidth',1);
patch(circle2(1,:), circle2(2,:), circle2(3,:),'r');

for i=1:n
    Laser_Pose(:,:,i)= INI*trotx(-pi/2)*troty(2*pi*i/n)*transl(0, 0, -radius);
end

Q= p560.ikine6s(Laser_Pose, 'run');
p560.plot(Q);