clc; close all; clear all
earth_image = imread('earth4.jpg');
%imshow(earth_image);
%input data
hours = 3600;
G = 6.6742e-20;
deg =pi/180;

%img 3 file size smaller easier to read
%global h1
% Earth data:
m1 = 5.974e24; %kg
R = 6378;  %km
m2 = 260;  %kg

G = 6.6742e-20; %si unit
mu = G*(m1 + m2);  %(km^3/s^2) gravitation parameneter
num_sat = 3; %%%%%CHANGE IF ADDING MORE %%%%


%satellite 1 data (Starlink 5257)
e1 =0.00036; %Eccentricity
RA1 = 309.4365; %deg (20.6291 hr) Right ascension
incl1 = 70.0013; %degrees
w1 =249.8617; %degrees   Argument of perigee
TA1 =0; %assume 0  True anomaly
a1 = R + 362.8; % semi major axis in km
h1 = sqrt(mu*a1*(1 - e1^2)); %(km^2/s) Angular momentum

%Enter satellite 2 data (Starlink-2227)
e2 =0.00010;
RA2= 1.05; %deg (0.070 hr)
incl2 = 53.0564; %degrees
w2 =59.0203	; %degrees
TA2 =0; %assume 0
a2 = R+554.2; % semi major axis in km
h2 = sqrt(mu*a2*(1 - e2^2)); %(km^2/s)

%Enter satellite 3 data (Starlink-5642)
e3 =0.00027;
RA3= 325.3845; %deg (21.6923  hr)
incl3 = 69.9991; %degrees
w3 =271.9722; %degrees
TA3 =0; %assume 0
a3 = R+576.8 ; % semi major axis in km
h3 = sqrt(mu*a3*(1 - e3^2)); %(km^2/s)


%     e4 =0.00049;
%     RA4= 333.2955; %deg (22.22 hr)
%     incl4 = 53.051	; %degrees
%     w4 =153.1184; %degrees
%     TA4 =0; %assume 0
%     a4 = R+383.2; % semi major axis in km
%     h4 = sqrt(mu*a4*(1 - e4^2)); %(km^2/s)
%     global mu, RA1; RA2; RA3;
%
%appending variables to make it easier to pass
coe ={[h1 ,e1 ,RA1*deg, incl1*deg, w1*deg ,TA1*deg],[h2 ,e2 ,RA2*deg ,incl2*deg, w2*deg, TA2*deg],[h3 ,e3 ,RA3*deg ,incl3*deg, w3*deg, TA3*deg]};
%,[h4 ,e4 ,RA4*deg, incl4*deg, w4*deg ,TA4*deg]};

%coe = reshape(coe, num_sat,[]); for unti 17 if i required (backup)

for j = 1 : num_sat
    [r, v] =sv_from_coe(coe(1,j), mu);

    r0 = [r(1) r(2) r(3)];
    v0 = [v(1) v(2) v(3)];
    t0 = 0;
    tf = 2 *hours; %source 90 minutes to 110 period
    y0 = [r0 v0]';
    [t,y] = ode45(@rates, [t0 tf], y0);

    %Output the results:

    for i = 1:length(t)
        r(i) = norm([y(i,1) y(i,2) y(i,3)]);
    end
    [rmax, imax] = max(r);
    [rmin, imin] = min(r);
    v_at_rmax = norm([y(imax,4) y(imax,5) y(imax,6)]);
    v_at_rmin = norm([y(imin,4) y(imin,5) y(imin,6)]);
    %...Output to the command window:
    fprintf('\n--------------------------------------------------------\n')
    fprintf('\n The initial position is [%g, %g, %g] (km).',...
        r0(1), r0(2), r0(3))
    fprintf('\n Magnitude = %g km\n', norm(r0))
    fprintf('\n The initial velocity is [%g, %g, %g] (km/s).',...
        v0(1), v0(2), v0(3))
    fprintf('\n Magnitude = %g km/s\n', norm(v0))
    fprintf('\n Initial time = %g h.\n Final time = %g h.\n',0,tf/hours)
    fprintf('\n The minimum altitude is %g km at time = %g h.',...
        rmin-R, t(imin)/hours)
    fprintf('\n The speed at that point is %g km/s.\n', v_at_rmin)
    fprintf('\n The maximum altitude is %g km at time = %g h.',...
        rmax-R, t(imax)/hours)
    fprintf('\n The speed at that point is %g km/s\n', v_at_rmax)
    fprintf('\n--------------------------------------------------------\n\n')
    %Plot the results:
    %Draw the planet
    [xx, yy, zz] = sphere(100);
    surf(R*xx, R*yy, R*-zz, R*earth_image)
    h = findobj('Type', 'surface');

    set(h, 'Cdata', earth_image, 'FaceColor' , 'texturemap', 'edgecolor', 'none')
    colormap(light_gray)
    clim([-R/100 R/100])

    %     shading interp %not needed anymore
    % Draw and label the X, Y and Z axes
    line([0 2*R], [0 0], [0 0]); text(2*R, 0, 0, 'X')
    line( [0 0], [0 2*R], [0 0]); text( 0, 2*R, 0, 'Y')
    line( [0 0], [0 0], [0 2*R]); text( 0, 0, 2*R, 'Z')
    %% Plotting orbit
    hold on
    plot3( y(:,1), y(:,2), y(:,3),'k')
%     line([0 r0(1)], [0 r0(2)], [0 r0(3)])
%     text( y(1,1), y(1,2), y(1,3), 'o') %inital starting position componenet
%     text( y(end,1), y(end,2), y(end,3), 'f') %final position component

    %% inital view
    view([1, 1,.4])

    %% graph properties
    grid on
    axis equal
    xlabel('km')
    ylabel('km')
    zlabel('km')

    hold on

end






