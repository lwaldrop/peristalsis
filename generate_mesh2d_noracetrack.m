clear all
close all

L = 1;                              % length of computational domain (m)
N = 512;                            % number of Cartesian grid meshwidths at the finest level of the AMR grid
dx = L/N;                           % Cartesian mesh width (m)
ds = L/(2*N);                       % space between boundary points in straight tube

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parameters for the racetrack

Let = 0.3;                          % Length of elastic tube (m)
Nend = 0;                          % Number of rigid points on each end of elastic section
Lt = Let+2*Nend*ds;                 % Length of straight section with three rigid points on each end
Llong = 0;                        % Lenght of long end of the tube

diameter = 0.1;                     % diameter of the tube
R2 = 0.1;                           % radius of inner wall
R1 = R2+diameter;                   % radius of outer wall

Nstraight = 2*ceil(Let/ds)          % number of points along each straight section
%Ncurve = 2*ceil(pi*R1/ds);          % number of points along each curved section
Nrace = 4*ceil(Llong/ds);         % number of points making up the racetrack part
%Nrace = Nstraight+2*Ncurve;         % number of points making up the racetrack part
%dtheta = pi/(Ncurve/2);             % angle increment for drawing curved edges

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
kappa_spring = 100;               % spring constant (Newton)
kappa_beam = 1.0;                 % beam stiffness constant (Newton m^2)
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
for i=0:ceil(Nstraight/2-1),
    ytop = centery-R2;
    xtop = -Lt/2+(i-1)*ds;
    plot(xtop,ytop,'r-')
    fprintf(vertex_fid, '%1.16e %1.16e\n', xtop, ytop);
end

%bottom part
for i=ceil(Nstraight/2):Nstraight-1,
    ybot = centery-R1;
    xbot = -Lt/2+(i-ceil(Nstraight/2)-1)*ds;
    plot(xbot,ybot,'r-')
    fprintf(vertex_fid, '%1.16e %1.16e\n', xbot, ybot);
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
        plot(x,y,'k.')
    fprintf(vertex_fid, '%1.16e %1.16e\n', x, y);
    end
end


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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NO race track part
% Write out the vertex information

vertex_fid = fopen([mesh_name 'race_' num2str(N) '.vertex'], 'w');
fprintf(vertex_fid, '%d\n', Nrace);

%right part of NO racetrack
for i=1:ceil(Nrace/4),
    ytop = centery-R2;
    xtop = Lt/2+i*ds;
    plot(xtop,ytop,'b-')
    fprintf(vertex_fid, '%1.16e %1.16e\n', xtop, ytop);
end

for i=1:ceil(Nrace/4),
    ybot = centery-R1;
    xbot = Lt/2+i*ds;
    plot(xbot,ybot,'b-')
    fprintf(vertex_fid, '%1.16e %1.16e\n', xbot, ybot);
end


%left part of NO racetrack
for i=1:ceil(Nrace/4),
    ytop = centery-R2;
    xtop = -Lt/2-i*ds;
    plot(xtop,ytop,'b-')
    fprintf(vertex_fid, '%1.16e %1.16e\n', xtop, ytop);
end

for i=1:ceil(Nrace/4),
    ybot = centery-R1;
    xbot = -Lt/2-i*ds;
    plot(xbot,ybot,'b-')
    fprintf(vertex_fid, '%1.16e %1.16e\n', xbot, ybot);
end
fclose(vertex_fid);




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
%target_fid = fopen([mesh_name 'tube_' num2str(N) '.target'], 'w');

%fprintf(target_fid, '%d\n', 4*Nend);

%for i = 0:Nend-1,
%    fprintf(target_fid, '%d %1.16e\n', i, kappa_target*ds/(ds^2));
%end

%for i = 0ceil(Nstraight/2):ceil(Nstraight/2)-1,
%    fprintf(target_fid, '%d %1.16e\n', i, kappa_target*ds/(ds^2));
%end

%for i = ceil(Nstraight/2):ceil(Nstraight/2)+Nend-1,
%    fprintf(target_fid, '%d %1.16e\n', i, kappa_target*ds/(ds^2));
%end

%for i = Nstraight-Nend:Nstraight-1,
%    fprintf(target_fid, '%d %1.16e\n', i, kappa_target*ds/(ds^2));
%end

%fclose(target_fid);


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
% Write out the target point information for the racetrack
target_fid = fopen([mesh_name 'race_' num2str(N) '.target'], 'w');

fprintf(target_fid, '%d\n', Nrace);

for i = 0:Nrace-1,
    fprintf(target_fid, '%d %1.16e\n', i, kappa_target*ds/(ds^2));
end

fclose(target_fid);
