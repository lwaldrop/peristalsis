close all
clear all

L = 1;                              % length of computational domain (m)
N = 512;                            % number of Cartesian grid meshwidths at the finest level of the AMR grid
dx = L/N;                           % Cartesian mesh width (m)
ds = L/(2*N);                       % space between boundary points in straight tube

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parameters for the racetrack

Let = 0.4;                          % Length of elastic tube (m)
Nend = 10;                          % Number of rigid points on each end of elastic section
Lt = Let+2*Nend*ds;                 % Length of straight section with three rigid points on each end
Llong = 0.05;                        % Lenght of long end of the tube

diameter = 0.1;                     % diameter of the tube
R2 = 0.1;                           % radius of inner wall
theta = 0.45;
R1 = R2+diameter;                   % radius of outer wall
Br1dia = diameter*nthroot(0.5,3)   % radius of first branch, according to Murray's law
x_addfactor = sqrt(-(diameter/2)^2+(Br1dia/cos(theta))^2);

Nstraight = 2*ceil(Lt/ds);          % number of points along each straight section
%Ncurve = 2*ceil(pi*R1/ds);          % number of points along each curved section
nseg = 10;                             %number of segments in the racetrack
Nrace = nseg*ceil(Llong/ds);         % number of points making up the racetrack part
%Nrace = Nstraight+2*Ncurve;         % number of points making up the racetrack part

mesh_name = 'heart_';               % structure name

centery = 0;                        % y-position of center of curved sections
centerx1 = -0.5*Lt;                 % x-position of center of left curved section
centerx2 = 0.5*Lt;                  % x-position of center of right curved section

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters for the pericardium
Dp = 2*diameter;                    %diameter of the pericardium
Nperi = 2*ceil((Dp-diameter)/ds);  % number of boundary points along the sides of the pericardium
Nperitot = Nperi + Nstraight;       % total number of pericardium points

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters for the actuator
La = 0.04;                          % length of the actuator section
NLa = ceil(La/ds);                  % number of points along each actuator
Ca = 0.25*Lt;                       % center of the actuator section
NCa = ceil(ceil(Nstraight/2)*Ca/Lt); % index of the center point
Na1 = NCa - ceil(NLa/2);            % index of the starting point with respect to the elastic section
Na2 = Na1+NLa-1;                    % index of the ending point with respect to the elastic section

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters for the prescribed peristalsis
Lap = 0.3;                          % length of the actuator section
NLap = ceil(Lap/ds);                  % number of points along each actuator
Cap = 0.5*Lt;                       % center of the actuator section
NCap = ceil(ceil(Nstraight/2)*Cap/Lt); % index of the center point
Na1p = NCap - ceil(NLap/2);            % index of the starting point with respect to the elastic section
Na2p = Na1p+NLap-1;                    % index of the ending point with respect to the elastic section


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters for markers
Nmarkersx = 11;                     %number of columns of markers
Nmarkersy = 11;                     %number of markers in each column
Nmarkers=Nmarkersx*Nmarkersy;       %total number of markers
dmx = Let/(Nmarkersx-1);            %space between markers in x-direction
dmy = diameter/(Nmarkersy-1);       %space between markers in y-direction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% material parameters
kappa_spring = 2.0e0;               % spring constant (Newton)
kappa_beam = 5.0e-1;                 % beam stiffness constant (Newton m^2)
kappa_target = kappa_spring;        % target point penalty spring constant (Newton)


%%%%%%%%%%%%%
% Plotting things 

figure(1)
hold on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make the elastic section of the tube
% Write out the vertex information

vertex_fid = fopen([mesh_name 'tube_' num2str(N) '.vertex'], 'w');
fprintf(vertex_fid, '%d\n', Nstraight);

%top part
for i=1:ceil(Nstraight/2),
    ytop = centery-R2;
    xtop = -Lt/2+(i-1)*ds;
    plot(xtop,ytop,'r-+')
    axis([-0.5 0.5 -0.5 0.5])
    fprintf(vertex_fid, '%1.16e %1.16e\n', xtop, ytop);
end

%bottom part
for i=ceil(Nstraight/2)+1:Nstraight,
    ybot = centery-R1;
    xbot = -Lt/2+(i-ceil(Nstraight/2)-1)*ds;
    plot(xbot,ybot,'r-+')
    fprintf(vertex_fid, '%1.16e %1.16e\n', xbot, ybot);
end
fclose(vertex_fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prescribed peristalsis part
% Write out the vertex information

%top part
vertex_fid = fopen(['pperi_top_' num2str(N) '.vertex'], 'w');
fprintf(vertex_fid, '%d\n', NLap);

for i=Na1p:Na2p,
    ytop = centery-R2;
    xtop = -Lt/2+i*ds;
    plot(xtop,ytop,'m*')
    fprintf(vertex_fid, '%1.16e %1.16e\n', xtop, ytop);
end
fclose(vertex_fid);

%bottom part
vertex_fid = fopen(['pperi_bot_' num2str(N) '.vertex'], 'w');
fprintf(vertex_fid, '%d\n', NLap);

for i=Na1p:Na2p,
    ybot = centery-R1;
    xbot = -Lt/2+i*ds;
    plot(xbot,ybot,'m*')
    fprintf(vertex_fid, '%1.16e %1.16e\n', xbot, ybot);
end
fclose(vertex_fid);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % actuator part
% % Write out the vertex information
% 
% %top part
% vertex_fid = fopen(['actuator_top_' num2str(N) '.vertex'], 'w');
% fprintf(vertex_fid, '%d\n', NLa);
% 
% for i=Na1:Na2,
%     ytop = centery-R2;
%     xtop = -Lt/2+i*ds;
%     plot(xtop,ytop,'k*')
%     fprintf(vertex_fid, '%1.16e %1.16e\n', xtop, ytop);
% end
% fclose(vertex_fid);
% 
% %bottom part
% vertex_fid = fopen(['actuator_bot_' num2str(N) '.vertex'], 'w');
% fprintf(vertex_fid, '%d\n', NLa);
% 
% for i=Na1:Na2,
%     ybot = centery-R1;
%     xbot = -Lt/2+i*ds;
%     plot(xbot,ybot,'k*')
%     fprintf(vertex_fid, '%1.16e %1.16e\n', xbot, ybot);
% end
% fclose(vertex_fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BRANCHED race track part
% Write out the vertex information

vertex_fid = fopen([mesh_name 'race_' num2str(N) '.vertex'], 'w');
fprintf(vertex_fid, '%d\n', 2*Nrace);

%right part of straight section
for i=1:ceil(Nrace/nseg),
    ytop = centery-R2;
    xtop = Lt/2+i*ds;
    plot(xtop,ytop,'b-')
    fprintf(vertex_fid, '%1.16e %1.16e\n', xtop, ytop);
end

for i=1:ceil(Nrace/nseg),
    ybot = centery-R1;
    xbot = Lt/2+i*ds;
    plot(xbot,ybot,'b-')
    fprintf(vertex_fid, '%1.16e %1.16e\n', xbot, ybot);
end

%right part of branched section

for i=1:ceil(2*(Nrace/nseg)),
    ytop = centery-R2+i*ds*sin(theta);
    xtop = Lt/2+Llong+i*ds;
    plot(xtop,ytop,'y-')
    fprintf(vertex_fid, '%1.16e %1.16e\n', xtop, ytop);
end

xtop1 = xtop;
ytop1 = ytop;

for i=1:ceil(2*(Nrace/nseg)),
    ybot = centery-R1-i*ds*sin(theta);
    xbot = Lt/2+Llong+i*ds;
    plot(xbot,ybot,'g-')
    fprintf(vertex_fid, '%1.16e %1.16e\n', xbot, ybot);
end

xbot2 = xbot;
ybot2 = ybot;

%right part of inner branched section
for i=1:ceil(2*(Nrace/nseg)),
    ytop = centery-R2-0.5*diameter+i*ds*sin(theta);
    xtop = Lt/2+Llong+x_addfactor+i*ds;
    plot(xtop,ytop,'b-')
    fprintf(vertex_fid, '%1.16e %1.16e\n', xtop, ytop);
end

xbot1 = xtop;
ybot1 = ytop;

for i=1:ceil(2*(Nrace/nseg)),
    ybot = centery-R2-0.5*diameter-i*ds*sin(theta);
    xbot = Lt/2+Llong+x_addfactor+i*ds;
    plot(xbot,ybot,'b-')
    fprintf(vertex_fid, '%1.16e %1.16e\n', xbot, ybot);
end

xtop2 = xbot;
ytop2 = ybot;


plot(xtop1,ytop1,'k*')
plot(xbot1,ybot1,'k*')
branch1topdiameter = cos(theta)*sqrt((xtop1-xbot1)^2+(ytop1-ybot1)^2)
plot(xtop2,ytop2,'m*')
plot(xbot2,ybot2,'m*')
branch1botdiameter = cos(theta)*sqrt((xtop2-xbot2)^2+(ytop2-ybot2)^2)

if round(branch1topdiameter) == round(Br1dia) && round(branch1botdiameter) == round(Br1dia)
    disp('This system branches according to Murrays law!')
else
    disp('WARNING!! WARNING!!')
    disp('This system does NOT branch according to Murrays law!')
end


%left part of straight section
for i=1:ceil(Nrace/nseg),
    ytop = centery-R2;
    xtop = -Lt/2-i*ds;
    plot(xtop,ytop,'b-')
    fprintf(vertex_fid, '%1.16e %1.16e\n', xtop, ytop);
end

for i=1:ceil(Nrace/nseg),
    ybot = centery-R1;
    xbot = -Lt/2-i*ds;
    plot(xbot,ybot,'b-')
    fprintf(vertex_fid, '%1.16e %1.16e\n', xbot, ybot);
end

%left part of outer branched section
for i=1:ceil(2*(Nrace/nseg)),
    ytop = centery-R2+i*ds*sin(theta);
    xtop = -Lt/2-Llong-i*ds;
    plot(xtop,ytop,'b-')
    fprintf(vertex_fid, '%1.16e %1.16e\n', xtop, ytop);
end


xtop1 = xtop;
ytop1 = ytop;

for i=1:ceil(2*(Nrace/nseg)),
    ybot = centery-R1-i*ds*sin(theta);
    xbot = -Lt/2-Llong-i*ds;
    plot(xbot,ybot,'b-')
    fprintf(vertex_fid, '%1.16e %1.16e\n', xbot, ybot);
end

xbot2 = xbot;
ybot2 = ybot;


%left part of inner branched section
for i=1:ceil(2*(Nrace/nseg)),
    ytop = centery-R2-0.5*diameter+i*ds*sin(theta);
    xtop = -Lt/2-Llong-x_addfactor-i*ds;
    plot(xtop,ytop,'b-')
    fprintf(vertex_fid, '%1.16e %1.16e\n', xtop, ytop);
end

xbot1 = xtop;
ybot1 = ytop;

for i=1:ceil(2*(Nrace/nseg)),
    ybot = centery-R2-0.5*diameter-i*ds*sin(theta);
    xbot = -Lt/2-Llong-x_addfactor-i*ds;
    plot(xbot,ybot,'b-')
    fprintf(vertex_fid, '%1.16e %1.16e\n', xbot, ybot);
end


xtop2 = xbot;
ytop2 = ybot;

plot(xtop1,ytop1,'k*')
plot(xbot1,ybot1,'k*')
branch1topdiameter = cos(theta)*sqrt((xtop1-xbot1)^2+(ytop1-ybot1)^2)
plot(xtop2,ytop2,'m*')
plot(xbot2,ybot2,'m*')
branch1botdiameter = cos(theta)*sqrt((xtop2-xbot2)^2+(ytop2-ybot2)^2)

if round(branch1topdiameter) == round(Br1dia) && round(branch1botdiameter) == round(Br1dia)
    disp('This system branches according to Murrays law!')
else
    disp('WARNING!! WARNING!!')
    disp('This system does NOT branch according to Murrays law!')
end


fclose(vertex_fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make markers as vertices with no material properties

vertex_fid = fopen(['markers_' num2str(N) '.vertex'], 'w');
fprintf(vertex_fid, '%d\n', Nmarkers);

%top part
for i=0:Nmarkersx-1,
    for j=0:Nmarkersy-1,
        y = centery-R2-j*dmy;
        x = -Let/2+i*dmx;
    fprintf(vertex_fid, '%1.16e %1.16e\n', x, y);
    plot(x,y)
    end
end
fclose(vertex_fid);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pericardium
% Write out the vertex information
% 
% vertex_fid = fopen([mesh_name 'peri_' num2str(N) '.vertex'], 'w');
% fprintf(vertex_fid, '%d\n', Nperitot);
% 
% % make the top and bottom of the pericardium
% for i=1:ceil(Nstraight/2),
%     ytop = centery-(R2-(Dp-diameter)/2);
%     xtop = -Lt/2+(i-1)*ds;
%     fprintf(vertex_fid, '%1.16e %1.16e\n', xtop, ytop);
% end
% 
% for i=ceil(Nstraight/2)+1:Nstraight,
%     ybot = centery-R1-(Dp-diameter)/2;
%     xbot = -Lt/2+(i-ceil(Nstraight/2)-1)*ds;
%     fprintf(vertex_fid, '%1.16e %1.16e\n', xbot, ybot);
% end
% 
% % make the four side pieces
% for i=Nstraight+1:Nstraight+ceil(Nperi/4),
%     y = centery-(R1+(Dp-diameter)/2)+(i-Nstraight-1)*ds;
%     x = -Lt/2;
%     fprintf(vertex_fid, '%1.16e %1.16e\n', x, y);
% end
% 
% for i=Nstraight+ceil(Nperi/4)+1:Nstraight+ceil(Nperi/2),
%     y = centery-R2+(i-Nstraight-ceil(Nperi/4)-1)*ds;
%     x = -Lt/2;
%     fprintf(vertex_fid, '%1.16e %1.16e\n', x, y);
% end
% 
% for i=Nstraight+ceil(Nperi/2)+1:Nstraight+ceil(3*Nperi/4),
%     y = centery-(R1+(Dp-diameter)/2)+(i-Nstraight-ceil(Nperi/2)-1)*ds;
%     x = Lt/2;
%     fprintf(vertex_fid, '%1.16e %1.16e\n', x, y);
% end
% 
% for i=Nstraight+ceil(3*Nperi/4)+1:Nperitot,
%     y = centery-R2+(i-Nstraight-ceil(3*Nperi/4)-1)*ds;
%     x = Lt/2;
%     fprintf(vertex_fid, '%1.16e %1.16e\n', x, y);
% end
% fclose(vertex_fid);


hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write out the spring information for the elastic section

spring_fid = fopen([mesh_name 'tube_' num2str(N) '.spring'], 'w');
fprintf(spring_fid, '%d\n', Nstraight-2);

%elastic part of tube
for i = 0:ceil(Nstraight/2)-2,
    fprintf(spring_fid, '%d %d %1.16e %1.16e\n', i, i+1, kappa_spring*ds/(ds^2), ds);
end

for i = ceil(Nstraight/2):Nstraight-2,
    fprintf(spring_fid, '%d %d %1.16e %1.16e\n', i, i+1, kappa_spring*ds/(ds^2), ds);
end

fclose(spring_fid);
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Write out the beam information for the elastic section

beam_fid = fopen([mesh_name 'tube_' num2str(N) '.beam'], 'w');
fprintf(beam_fid, '%d\n', Nstraight-4);

%elastic part of tube
for i = 0:ceil(Nstraight/2)-3,
    fprintf(beam_fid, '%d %d %d %1.16e\n', i, i+1, i+2, kappa_beam*ds/(ds^4));
end

for i = ceil(Nstraight/2):Nstraight-3,
    fprintf(beam_fid, '%d %d %d %1.16e\n', i, i+1, i+2, kappa_beam*ds/(ds^4));
end
fclose(beam_fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Write out the target point information for the ends of the elastic tube
target_fid = fopen([mesh_name 'tube_' num2str(N) '.target'], 'w');

fprintf(target_fid, '%d\n', 4*Nend);

for i = 0:Nend-1,
    fprintf(target_fid, '%d %1.16e\n', i, kappa_target*ds/(ds^2));
end

for i = ceil(Nstraight/2)-Nend:ceil(Nstraight/2)-1,
    fprintf(target_fid, '%d %1.16e\n', i, kappa_target*ds/(ds^2));
end

for i = ceil(Nstraight/2):ceil(Nstraight/2)+Nend-1,
    fprintf(target_fid, '%d %1.16e\n', i, kappa_target*ds/(ds^2));
end

for i = Nstraight-Nend:Nstraight-1,
    fprintf(target_fid, '%d %1.16e\n', i, kappa_target*ds/(ds^2));
end

fclose(target_fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write out the target point information for the actuator

%top actuator
target_fid = fopen(['actuator_top_' num2str(N) '.target'], 'w');
fprintf(target_fid, '%d\n', NLa);

for i = 0:NLa-1,
    fprintf(target_fid, '%d %1.16e\n', i, kappa_target*ds/(ds^2));
end

fclose(target_fid);

%bottom actuator
target_fid = fopen(['actuator_bot_' num2str(N) '.target'], 'w');
fprintf(target_fid, '%d\n', NLa);

for i = 0:NLa-1,
    fprintf(target_fid, '%d %1.16e\n', i, kappa_target*ds/(ds^2));
end

fclose(target_fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write out the target point information for the racetrack
target_fid = fopen([mesh_name 'race_' num2str(N) '.target'], 'w');

fprintf(target_fid, '%d\n', 2*Nrace);

for i = 0:2*Nrace-1,
    fprintf(target_fid, '%d %1.16e\n', i, kappa_target*ds/(ds^2));
end

fclose(target_fid);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write out the target point information for the actuator

%top prescribed peristalsis
target_fid = fopen(['pperi_top_' num2str(N) '.target'], 'w');
fprintf(target_fid, '%d\n', NLap);

for i = 0:NLap-1,
    fprintf(target_fid, '%d %1.16e\n', i, kappa_target*ds/(ds^2));
end

fclose(target_fid);

%bottom prescribed peristalsis
target_fid = fopen(['pperi_bot_' num2str(N) '.target'], 'w');
fprintf(target_fid, '%d\n', NLap);

for i = 0:NLap-1,
    fprintf(target_fid, '%d %1.16e\n', i, kappa_target*ds/(ds^2));
end

fclose(target_fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write out the target point information for the pericardium
% target_fid = fopen([mesh_name 'peri_' num2str(N) '.target'], 'w');
% 
% fprintf(target_fid, '%d\n', Nperitot);
% 
% for i = 0:Nperitot-1,
%     fprintf(target_fid, '%d %1.16e\n', i, kappa_target*ds/(ds^2));
% end
% 
% fclose(target_fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
