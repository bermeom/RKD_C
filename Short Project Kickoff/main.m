%% Lab1: ROBOTICS , KINEMATICS, DYNAMICS AND CONTROL
% Miguel Angle Bermeo Ayerbe & Guillermo Lopez Diest    
close all;clear all;clc;

%% The Robot task
% Drilling holes in rear and front faces of an iron piece with shape of 
% half torus. Milling a groove-shape around the outer radius of the torus 
% and finally weld the iron strip joints which conform a spiral with a 
% half-torus shape The main dimensions of the torus are: 0.3m of tube 
% diameter and 0.95 m outer radius. A specialized machine conform the 
% tube-spiral in a torus shape. See a video of a spiral tubeformer 
% machinery: 
%%
% <https://www.youtube.com/watch?v=1rvbirDwvRI>
%
%

% This is some data you can use
load('Data_groove_weld_fv_torus.mat');
% Rotate the point of view to better understand the task to be done
open('Torus.fig') 

%% Task 2. Drilling and milling 8 shapes
% Numers of groove 8 in the outer diameter The groove has the following shape in milimeteres.
%
figure
plot3(Groove(1,:),Groove(2,:),Groove(3,:),'r') % plotting the Groove
title ('Groove shape')
xlabel('x');
ylabel('y');
zlabel('z');
axis equal;
%% Task 3 Welding points
% The Torus part is conformed from a iron strip that generate a spiral-like 
% tube. This assembling must be weld to fixe the structure and gain in 
% rigidity. The manipulation of the machine outputs a Torus of 8 turns and 
% dimensions: 0.3m of tube diameter and 0.95 m outer radius It is planned 
% to first weld every 117.8 mm. See the solution for better understanding.
%
%
figure
plot3(Weld_points(1,:),Weld_points(2,:),Weld_points(3,:),'g','LineWidth',2)
hold on
scatter3(Weld_points(1,:),Weld_points(2,:),Weld_points(3,:),'b')%;,'fillet');
xlabel('x');
ylabel('y');
zlabel('z');
axis 'equal';
grid on;