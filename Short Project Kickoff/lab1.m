%% Lab1: ROBOTICS , KINEMATICS, DYNAMICS AND CONTROL
% Miguel Angle Bermeo Ayerbe
% install startup_rvc
% startup_rvc
close all;clear all;clc; 

%% Algorithm
clear all; close all;clc; 
d=0.3;
n=8;
m=8*n;
r = 0.95;
r1 = 0.95-d/2;
r2 = d/2+0.001;
scale = 0.2;
hr = 2.5;
Ptos1 =[];
for i=0:(m)
    Ptos1(:,:,i+1)=troty((-pi*i)/m)*transl(r1,0,0)*trotz((2*pi*i)/n)*transl(r2,0,0);
end
spiral=transl(Ptos1)';
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
Ptos2 = [];
for i=0:(m1-1)
    Ptos2(:,:,i*2+1)=troty((-pi*(i+1)+pi/2)/m1)*transl(r1,0,0)*trotz(-pi/2)*transl(r2,0,0);
    Ptos2(:,:,i*2+2)=troty((-pi*(i+1)+pi/2)/m1)*transl(r1,0,0)*trotz(pi/2)*transl(r2,0,0);
end
holes1=transl(Ptos2)';

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
Ptos3=[];
for i=1:(m2)
        for (j=1:n2)
            Ptos3(:,:,(i-1)*n2+j)=troty((-pi*(i))/(m2)+pi/24+(pi*j)/100)*transl(r,0,0);%*trotz(2*pi)*transl(r2,0,0);
        end
end
holes2=transl(Ptos3)';

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

w=2*0.6;
h=2*1.2;
th1=pi*(20/180); 
th2=pi*(30/180);
tr1=0.75;
T1 = trotx(th2)*troty(th1);
T4 = transl(h/2,w/2,0);
T3 = T4*transl(r1,0,0);
table = [0 h h 0;...
         0 0 w w;...
         0 0 0 0]; 
temp = [table;ones(1,size(table,2))];
temp = T1*temp;
T2 = transl(-h/2,0,tr1+abs(min(temp(3,:))));
temp = T2*temp;
table = temp(1:3,:) 
figure
fill3(table(1,:),table(2,:),table(3,:),'r');
xlabel('x');
ylabel('y');
zlabel('z');
axis 'equal';
alpha 0.3;
hold on;


theta = pi*(30/180);%% 30
fv=stlread('Torus.stl');% fv is a struct with faces and vertices
fh = [fv.vertices'; ones(1,size(fv.vertices,1))];
fh=troty(pi)*fh;
fh(1:3,:) = (r).*(fh(1:3,:)./(max(max(fh(1:3,:)))));
fh = T1*T3*fh;
fh = T2*fh;
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
% view([-135 35]);



% Puma 560
% hold on

% L(1) = Link([0, 0, 0, pi/2],'modified');
% L(2) = Link([0, 0, 4.318*scale, 0],'modified');
% L(3) = Link([0, 1.5005*scale, 0.203*scale, -pi/2],'modified');
% L(4) = Link([0, 4.331*scale, 0, pi/2],'modified');
% L(5) = Link([0, 0, 0, -pi/2],'modified');
% L(6) = Link([0, 0, 0, 0],'modified');
% p560m_RKDC = SerialLink(L,'base',transl(0,0,hr),'name','Puma 560');
% p560m_RKDC
% qz = [0 0 0 0 0 0]
% p560m_RKDC.plot(qz,'floorlevel',0)


Ptos2_ = T1*T4.*Ptos2;
Ptos2_ = T2.*Ptos2_;

temp = [holes2;ones(1,size(holes2,2))];
temp= T1*T4*temp;
temp = T2*temp;
holes2_ = temp(1:3,:);
scatter3(holes2_(1,:),holes2_(2,:),holes2_(3,:),'m','filled');
xlabel('x');
ylabel('y');
zlabel('z');
grid on;
axis equal;


% Plot a sketch of the environment of the robot
mdl_puma560
p560.base= transl(0,0,hr);

p560.a
p560.plot(qz,'floorlevel',-1);
hold on;

%%
%%
Q= p560.ikine6s(Ptos2_, 'run');
p560.plot(Q);
%% Animating the Puma manipulator performing a task
% *Initialize the workspace*:
%%
% We clear the workspace and set the value of some constants of the problem
% clc;clear;close all;
radius=0.20;
n=200;

% Plot a sketch of the environment of the robot
mdl_puma560
% figure(6)
% p560.alhpa= p560.alpha*2;
p560.base= transl(0,0,hr);
% p560.d*scale;

p560.a
p560.plot(qz,'floorlevel',-1);
hold on;
%%
circle1 = circle([-0.25 0.25 0], radius);
circle2 = circle([-0.25 0.25 0], 0.5*radius);
plot3(circle1(1,:), circle1(2,:), circle1(3,:),'g','LineWidth',1);
patch(circle2(1,:), circle2(2,:), circle2(3,:),'r');
INI = transl(-0.25, 0.25,0); %center of the part
for i=1:n
    Laser_Pose(:,:,i)= INI*trotx(-pi/2)*troty(2*pi*i/n)*transl(0, 0, -radius);
end

%%
Q= p560.ikine6s(Laser_Pose, 'run');
p560.plot(Q);

%%
Q= p560m_RKDC.ikine6s(Laser_Pose, 'run');
p560m_RKDC.plot(Q);
