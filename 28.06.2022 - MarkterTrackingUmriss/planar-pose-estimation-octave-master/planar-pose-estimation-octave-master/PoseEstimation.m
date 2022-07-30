clear; 

%load packages
pkg load image;
pkg load geometry;



%define Camera matrix
f = 917;
px = 500;
py = 375;
K = [f 0 px; 0 f py; 0 0 1];



%load image
img = imread('marker_img_k.jpg');
figure(1);
imshow(img);
hold on;
title("marker image");



%compute image threshold
level = graythresh (img, 'moments');



%convert image to binary image
b_img = im2bw (img, level);
%invert binary image
b_img = ~b_img;
figure(2);
imshow(b_img);
hold on;
title("binary marker image");


%find boundaries
bs = bwboundaries (b_img, 8, 'noholes');
%discard small boundaries
i = 0;
for k = 1:numel (bs)
  if size(bs{k}, 1) > size(img, 2) / 5
    i = i+1;
  endif
endfor
bounds = cell(i,1);
j = 1;
for k = 1:numel (bs)
  if size(bs{k}, 1) > size(img, 2) / 5
    bounds{j} = bs{k};
    j = j+1;
  endif
endfor
%show boundaries
figure(4);
imshow(b_img);
hold on;
for k = 1:numel (bounds)
  plot (bounds {k} (:, 2), bounds {k} (:, 1), 'r', 'linewidth', 2);
endfor
title("boundaries");



%simplify polygons
simp_bounds = cell(i,1);
num4p = 0;
for k = 1:numel (bounds)
  simp_bounds{k} = simplifyPolyline (bounds{k}, 'tol', 2);
  %figure();
  %plot (simp_bounds {k} (:, 2), simp_bounds {k} (:, 1), 'r', 'linewidth', 2);
  hold on;
  title("simplified boundaries");
  if size(simp_bounds{k},1) == 5
    num4p = num4p+1;
  endif
endfor



%only polygons with 4 points are relevant (test for 5 because start/end point appears twice)
bounds_4p = cell(num4p,1);
j = 1;
for k = 1:numel(simp_bounds)
  if size(simp_bounds{k},1) ==5
    bounds_4p{j} = simp_bounds{k};
    j = j+1;
  endif
endfor



%polygons have to be convex (number of points must remain 5)
marker_candidate =  zeros(5,2);
for k = 1:numel (bounds_4p)
  conv = convexHull(bounds_4p{k});
  if size(conv, 1) == 5
    marker_candidate = bounds_4p{k};
    break;
  endif
endfor
figure(5);
plot (marker_candidate (:, 2), marker_candidate (:, 1), 'g', 'linewidth', 2);
hold on;
title("marker candidate");



%define marker model for projecting with K into first image using z=2*f
depth = 2*f;
width = 50
height = 50
XWFTL = [-width/2;-height/2;depth];
XWFTR = [width/2;-height/2;depth];
XWFBR = [width/2;height/2;depth];
XWFBL = [-width/2;height/2;depth];
XWF = [XWFTL XWFTR XWFBR XWFBL];



%project the 3D marker position into the initial image
x0h = K*XWF;
x0 = points2DDivideW(x0h)



%define marker model again in z=0 plane
XW0TL = [-width/2;-height/2; 0]; XW0TLH = [XW0TL ;1];
XW0TR = [width/2;-height/2; 0]; XW0TRH = [XW0TR ;1];
XW0BR = [width/2;height/2; 0]; XW0BRH = [XW0BR ;1];
XW0BL = [-width/2;height/2; 0]; XW0BLH = [XW0BL ;1];
XW = [XW0TL XW0TR XW0BR XW0BL];
XWH = [XW0TLH XW0TRH XW0BRH XW0BLH];



%compute homography between image marker coordinates and 3D marker model in z=0 plane
%only need x and y coordinates of marker model
XWH_no_z =  XWH; 
XWH_no_z(3, :) = [];
HW0 = getHomography(XWH_no_z, x0h);



%compute homography between image marker coordinates and coordinates from real image
%reshape representation of coordinates
x1 = ones(3,4);
x1(1,:) = marker_candidate(1:4,2)';
x1(2,:) = marker_candidate(1:4,1)';
H01 = getHomography(x0, x1);



%compute HW1
HW1 = H01 * HW0;




%extract R and T from the homography
KinvH = inv(K) * HW1;
% normalize
l = sqrt(norm(KinvH(:,1)) * norm(KinvH(:,1)))
R1 = KinvH(:,1)/l;
R2 = KinvH(:,2)/l;
T = KinvH(:,3)/l;
c = R1+R2;
p = cross(R1, R2);
d = cross(c, p);

R1b = 1/sqrt(2) * (c/norm(c) + d/norm(d));
R2b = 1/sqrt(2) * (c/norm(c) - d/norm(d));
R3 = cross(R1b, R2b);
RT2 = [R1b, R2b, R3, T]
RT = [KinvH(:,1) KinvH(:,2) cross(KinvH(:,1), KinvH(:,2)) KinvH(:,3)];



%create final projection matrix
%P = K*RT;
P = K*RT2;


%test if projection matrix works
x1
x1_projected = points2DDivideW(P * XWH)



%marker origin
OW = [0,0,0]; OWH = [OW, 1];
o1h = P * OWH'; o1 = points2DDivideW(o1h)
%define coordinate cross
CCXW = [width/3, 0, 0]; CCXWH = [CCXW, 1];
CCYW = [0, width/3, 0]; CCYWH = [CCYW, 1];
CCZW = [0, 0, -width/3]; CCZWH = [CCZW, 1];
%project coordinate cross
ccx1h = P * CCXWH'; ccx1 = points2DDivideW(ccx1h);
ccy1h = P * CCYWH'; ccy1 = points2DDivideW(ccy1h);
ccz1h = P * CCZWH'; ccz1 = points2DDivideW(ccz1h);



%show result
figure(6);

%3D marker model
%subplot(1,3,1)
subplot(1,2,1)
scatter3(XW(1,:), XW(2,:), XW(3,:))
hold on
plot3([XW(1,:) XW(1,1)] , [XW(2,:) XW(2,1)], [XW(3,:) XW(3,1)])
%coordinate cross
plot3([OW(1) CCXW(1)], [OW(2) CCXW(2)], [OW(3) CCXW(3)], "r", "linewidth", 1); % X
plot3([OW(1) CCYW(1)], [OW(2) CCYW(2)], [OW(3) CCYW(3)], "g", "linewidth", 1); % Y
plot3([OW(1) CCZW(1)], [OW(2) CCZW(2)], [OW(3) CCZW(3)], "b", "linewidth", 1); % Z
axis([-width, width, -height, height, -width, width], "square")
title("Initial marker position 3D")
hold off

%marker candidate
%subplot(1,3,2)
%scatter(x1(1,:), x1(2,:))
%hold on

%plot([x1(1,:) x1(1,1)], [x1(2,:) x1(2,1)], "r")
%coordinate cross
%plot([o1(1) ccx1(1)], [o1(2) ccx1(2)], "r", "linewidth", 1); % X
%plot([o1(1) ccy1(1)], [o1(2) ccy1(2)], "g", "linewidth", 1); % Y
%plot([o1(1) ccz1(1)], [o1(2) ccz1(2)], "b", "linewidth", 1); % Z
%axis([0, 2*px, 0, 2*py], "square")
%title("Marker candidate position")


%real image 
%subplot(1,3,3)
subplot(1,2,2)
imshow(img)
hold on
plot([x1(1,:) x1(1,1)], [x1(2,:) x1(2,1)], "r")
%coordinate cross
plot([o1(1) ccx1(1)], [o1(2) ccx1(2)], "r", "linewidth", 2); % X
plot([o1(1) ccy1(1)], [o1(2) ccy1(2)], "g", "linewidth", 2); % Y
plot([o1(1) ccz1(1)], [o1(2) ccz1(2)], "b", "linewidth", 2); % /
axis([0, 2*px, 0, 2*py], "square")
title("Marker candidate position")





  
  