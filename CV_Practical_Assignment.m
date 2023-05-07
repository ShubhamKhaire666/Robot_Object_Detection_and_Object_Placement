clc
clear all
close all


%% Calibration

p02=[907 426;849 429;920 463;858 465;1200 410;1249 407;1218 445;1184 453;523 505;517 565;780 334;913 508;1270 441;1204 486;930 326]';
p03=[-6.75 -21.55 5;-11.5 -21.55 5;-6.55 -26.35 5;-11.55 -26.15 5;20.35 -21.55 5;25.45 -21.55 5;20.15 -26.65 5;20.25 -21.75 0;-35 -30.85 7.5;-35.1 -36.55 7.5;-15.2 0.3 0;-6.65 -26.25 0;25.35 -26.55 5;20.45 -26.45 0;-0.1 0 0]';

points2d=p02;points3d=p03;
points3d(4,:)=1;
x_mean = mean(points2d(1,:));
y_mean = mean(points2d(2,:));

d_mean = sum(sqrt(((points2d(1,:)-x_mean).^2 + (points2d(2,:)-y_mean).^2)))/length(points2d);

T = [sqrt(2)/d_mean 0 -sqrt(2)*x_mean/d_mean;
    0 sqrt(2)/d_mean -sqrt(2)*y_mean/d_mean;
    0 0 1];

x_mean = mean(points3d(1,:));
y_mean = mean(points3d(2,:));
z_mean = mean(points3d(3,:));

D_mean = sum(sqrt(((points3d(1,:)-x_mean).^2 + (points3d(2,:)-y_mean).^2) + (points3d(3,:)-z_mean).^2))/length(points2d);

U = [sqrt(3)/D_mean 0 0 -sqrt(3)*x_mean/D_mean;
    0 sqrt(3)/D_mean 0 -sqrt(3)*y_mean/D_mean;
    0 0 sqrt(3)/D_mean -sqrt(3)*z_mean/D_mean;
    0 0 0 1];

norm_p3 = U * points3d;

points2d(3,:) = 1;
norm_p2 = T * points2d;

norm_M = calibrate(norm_p2,norm_p3);

M = inv(T) * norm_M'*U;

reproject_err=reprojection_error(M,p03,p02)


%%  Cube and Robot Detection

cube_list=["Red"]
for itr=1:size(cube_list,2)
% img_direction=input('Please enter image direction: ')
img=imread("calibration2.png");
img = imgaussfilt(img,1);
Z_cube=5;
Z_robot=7.5;
Z_circle=0;
robot_l=1.2;


red=img(:,:,1)-max(img(:,:,2),img(:,:,3))>35 & img(:,:,2)<30 ;
red=uint8(red)*255;
green=img(:,:,2)-max(img(:,:,1),img(:,:,3))>20;
green=uint8(green)*255;
blue=img(:,:,3)-max(img(:,:,1),img(:,:,2))>6 & img(:,:,1)<70;
blue=uint8(blue)*255;


% Remove any small objects
bw_r = bwareaopen(red, 110);
%         figure(2)
%         imshow(bw_r)
%     % %
bw_g = bwareaopen(green, 50);
%     figure(2)
%     imshow(bw_g)
%
bw_b = bwareaopen(blue,100);
%     figure(3)
%     imshow(bw_b)

stats_r=regionprops('table',bw_r,'Centroid','Area','MajorAxisLength','MinorAxisLength','Circularity',...
    'Eccentricity');
stats_g=regionprops('table',bw_g,'Centroid','Area','MajorAxisLength','MinorAxisLength','Circularity',...
    'Eccentricity');
stats_b=regionprops('table',bw_b,'Centroid','Area','MajorAxisLength','MinorAxisLength','Circularity',...
    'Eccentricity');

blueRange = [0.50 0.6 0.6 0.75 0.38 0.65];
redRange = [0.88 0.95 0.45 93 0.42 0.6];

% Convert the image from RGB to HSV color space
hsvImage = rgb2hsv(img);

% Extract the hue, saturation, and value channels
hue = hsvImage(:,:,1);
saturation = hsvImage(:,:,2);
value = hsvImage(:,:,3);

% Create binary masks of the blue and red circles
blueMask = (hue >= blueRange(1) & hue <= blueRange(2)) ...
    & (saturation >= blueRange(3) & saturation <= blueRange(4)) ...
    & (value >= blueRange(5) & value <= blueRange(6));
blueMask = bwareaopen(blueMask, 40);


redMask = (hue >= redRange(1) & hue <= redRange(2)) ...
    & (saturation >= redRange(3) & saturation <= redRange(4)) ...
    & (value >= redRange(5) & value <= redRange(6));
redMask = bwareaopen(redMask, 30);


% Find the centroids of the blue and red circles
blueProps = regionprops('table',blueMask, 'Centroid','Circularity');
[b,index]=sort(blueProps.Circularity);
blueCentroids=blueProps.Centroid(index(end),:);

redProps = regionprops('table',redMask, 'Centroid','Circularity');
[b,index]=sort(redProps.Circularity);
redCentroids=redProps.Centroid(index(end),:);

% Calculate the Euclidean distance between the centroids of the blue and red circles
distance = pdist2(blueCentroids, redCentroids);
distance2=reshape(distance,[size(blueCentroids,1) size(redCentroids,1)]);

% Find the indices of the blue and red centroids that have a distance between 70 to 100
indices = find(distance >= 30 & distance <= 140);
if isempty(indices)
    blueMask=bwareaopen(blueMask,70);
    blueProps = regionprops('table',blueMask, 'Centroid','Circularity');
    [b,index]=sort(blueProps.Circularity);
    blueCentroids=blueProps.Centroid(index(end),:);
    distance = pdist2(blueCentroids, redCentroids);
    distance2=reshape(distance,[size(blueCentroids,1) size(redCentroids,1)]);
    indices = find(distance >= 30 & distance <= 140);
end
% Calculate the Euclidean distance between the centroids of the blue and red circles
[blueIndex,redIndex]=find(distance2==distance(indices(1)));


if(height(blueCentroids) ~= 1)
    blueCentroidsInRange = blueCentroids(blueIndex,:);
else
    blueCentroidsInRange = blueCentroids;
end
if(height(redCentroids) ~= 1)
    redCentroidsInRange = redCentroids(redIndex,:);
else
    redCentroidsInRange = redCentroids;
end
% Arrow Detection
%     [robot_b,robot_r]=arrow(img)

%    Blue and red points of robot
robot_r = redCentroids;
robot_b = blueCentroids;

% Display the original image with the blue and red circle centroids and the distance between them
figure(1)
imshow(img)
hold on

% Add a red circle around each red circle centroid that is in the range
radius = 20;
if robot_r == redCentroidsInRange
    viscircles(robot_r, radius, 'Color', 'r');
end

% Add a blue circle around each blue circle centroid that is in the range
radius = 20;
if intersect(robot_b,blueCentroidsInRange)
    viscircles(robot_b, radius, 'Color', 'b');
end

% Find centroids of cubes and circles
[b,ind_r]=sort(stats_r.Area);
if isempty(ind_r)
    disp("Can not find the red cube")
    %         break
else
    centroid_2Dr=stats_r.Centroid(ind_r(end),:);
end
if length(ind_r)==1
    circle_2Dr=stats_r.Centroid(ind_r(end),:);
else
    circle_2Dr=stats_r.Centroid(ind_r(end-1),:);
end

[b,ind_g]=sort(stats_g.Area);
if isempty(ind_g)
    disp("Can not find the green cube")
    %         break
else
    centroid_2Dg=stats_g.Centroid(ind_g(end),:);
end
if length(ind_g)==1
    circle_2Dg=stats_g.Centroid(ind_r(end),:);
else
    circle_2Dg=stats_g.Centroid(ind_g(end-1),:);
end

[b,ind_b]=sort(stats_b.Area);
if isempty(ind_b)
    disp("Can not find the blue cube")
    %         break
else
    centroid_2Db=stats_b.Centroid(ind_b(end),:);
end
if length(ind_b)==1
    circle_2Db=stats_b.Centroid(ind_b(end),:);
else
    circle_2Db=stats_b.Centroid(ind_b(end-1),:);
end

% p2d=[circle_2Dr;circle_2Dg;circle_2Db;centroid_2Dr;centroid_2Dg;centroid_2Db;robot_r;robot_b]';
% p3d=[-15 0 0;0 0 0;15 0 0;-9 -24 2.5;22.5 -24 2.5;9 42 2.5;-35 -30.85 7.5;-35.1 -36.55 7.5]';

% Find 3D points of robot, cubes , and circles


A=[M(:,1) M(:,2) -[robot_r 1]' (Z_robot*M(:,3)+M(:,4))];
[U,S,V]=svd(A);
ss=V(:,end);
robot_r3d=ss(1:2)./ss(end);
A=[M(:,1) M(:,2) -[robot_b 1]' (Z_robot*M(:,3)+M(:,4))];
[U,S,V]=svd(A);
ss=V(:,end);
robot_b3d=ss(1:2)./ss(end);


A=[M(:,1) M(:,2) -[centroid_2Dr 1]' ((Z_cube/2)*M(:,3)+M(:,4))];
[U,S,V]=svd(A);
ss=V(:,end);
cube_r=ss(1:2)./ss(end);

A=[M(:,1) M(:,2) -[circle_2Dr 1]' ((Z_circle/2)*M(:,3)+M(:,4))];
[U,S,V]=svd(A);
ss=V(:,end);
circle_r=ss(1:2)./ss(end);

A=[M(:,1) M(:,2) -[centroid_2Dg 1]' ((Z_cube/2)*M(:,3)+M(:,4))];
[U,S,V]=svd(A);
ss=V(:,end);
cube_g=ss(1:2)./ss(end);

A=[M(:,1) M(:,2) -[circle_2Dg 1]' ((Z_circle/2)*M(:,3)+M(:,4))];
[U,S,V]=svd(A);
ss=V(:,end);
circle_g=ss(1:2)./ss(end);

A=[M(:,1) M(:,2) -[centroid_2Db 1]' ((Z_cube/2)*M(:,3)+M(:,4))];
[U,S,V]=svd(A);
ss=V(:,end);
cube_b=ss(1:2)./ss(end);

A=[M(:,1) M(:,2) -[circle_2Db 1]' ((Z_circle/2)*M(:,3)+M(:,4))];
[U,S,V]=svd(A);
ss=V(:,end);
circle_b=ss(1:2)./ss(end);


%% Robot Movement Function
% for itr=1:size(cube_list,2)
    % Red Cube
    if strcmp(cube_list(itr),'Red')
        v1=[robot_b3d(1)-robot_r3d(1),robot_b3d(2)-robot_r3d(2)];
        v2=[cube_r(1)-robot_b3d(1),cube_r(2)-robot_b3d(2)];
        dot_product=dot(v1,v2);
        v1_norm=norm(v1);
        v2_norm=norm(v2);
        theta_r=acosd(dot_product/(v1_norm*v2_norm));
        cross_p=cross([v1,0], [v2,0]);
        if cross_p(3) > 0
            theta_r = -theta_r;
        end

        v1=[cube_r(1)-robot_b3d(1),cube_r(2)-robot_b3d(2)];
        v2=[circle_r(1)-cube_r(1),circle_r(2)-cube_r(2)];
        dot_product=dot(v1,v2);
        v1_norm=norm(v1);
        v2_norm=norm(v2);
        theta_cr=acosd(dot_product/(v1_norm*v2_norm));
        cross_p=cross([v1,0], [v2,0]);
        if cross_p(3) > 0
            theta_cr = -theta_cr;
        end

        d_tocube=sqrt((cube_r(1)-robot_b3d(1))^2+(cube_r(2)-robot_b3d(2))^2);
        d_tocircle=sqrt((cube_r(1)-circle_r(1))^2+(cube_r(2)-circle_r(2))^2);
        fprintf('turn(%f) ; go(%f) ; grab() ; turn(%f) ; go(%f) ; let_go() ; go(-10)\n',theta_r,d_tocube,theta_cr,d_tocircle-12)

%         d_10=(circle_r(1)-cube_r(1))^2+(circle_r(2)-cube_r(2))^2;
%         t_b=sqrt((10+12.5)^2/d_10);
%         t_r=sqrt((10+11.4+6)^2/d_10);
%         robot_b3d=[circle_r(1), circle_r(2)] + t_b * [cube_r(1)-circle_r(1), cube_r(2)-circle_r(2)];
%         robot_r3d=[circle_r(1), circle_r(2)] + t_r * [cube_r(1)-circle_r(1), cube_r(2)-circle_r(2)];
%         newp_b=M*[robot_b3d';7.5;1];
%         robot_b=newp_b./newp_b(3);
%         newp_r=M*[robot_r3d';7.5;1];
%         robot_r=newp_r./newp_r(3);

        % Green Cube
    elseif strcmp(cube_list(itr),'Green')
        v1=[robot_b3d(1)-robot_r3d(1),robot_b3d(2)-robot_r3d(2)];
        v2=[cube_g(1)-robot_b3d(1),cube_g(2)-robot_b3d(2)];
        dot_product=dot(v1,v2);
        v1_norm=norm(v1);
        v2_norm=norm(v2);
        theta_g=acosd(dot_product/(v1_norm*v2_norm));
        cross_p=cross([v1,0], [v2,0]);
        if cross_p(3) > 0
            theta_g = -theta_g;
        end

        v1=[cube_g(1)-robot_b3d(1),cube_g(2)-robot_b3d(2)];
        v2=[circle_g(1)-cube_g(1),circle_g(2)-cube_g(2)];
        dot_product=dot(v1,v2);
        v1_norm=norm(v1);
        v2_norm=norm(v2);
        theta_cg=acosd(dot_product/(v1_norm*v2_norm));
        cross_p=cross([v1,0], [v2,0]);
        if cross_p(3) > 0
            theta_cg = -theta_cg;
        end
        d_tocube=sqrt((cube_g(1)-robot_b3d(1))^2+(cube_g(2)-robot_b3d(2))^2);
        d_tocircle=sqrt((cube_g(1)-circle_g(1))^2+(cube_g(2)-circle_g(2))^2);
        fprintf('turn(%f) ; go(%f) ; grab() ; turn(%f) ; go(%f) ; let_go() ; go(-10)\n',theta_g,d_tocube,theta_cg,d_tocircle-12)

%         d_10=(circle_g(1)-cube_g(1))^2+(circle_g(2)-cube_g(2))^2;
%         t_b=sqrt((10+13.5)^2/d_10);
%         t_r=sqrt((10+11.4+6)^2/d_10);
%         robot_b3d=[circle_g(1), circle_g(2)] + t_b * [cube_g(1)-circle_g(1), cube_g(2)-circle_g(2)];
%         robot_r3d=[circle_g(1), circle_g(2)] + t_r * [cube_g(1)-circle_g(1), cube_g(2)-circle_g(2)];
%         newp_b=M*[robot_b3d';7.5;1];
%         robot_b=newp_b./newp_b(3);
%         newp_r=M*[robot_r3d';7.5;1];
%         robot_r=newp_r./newp_r(3);

        % Blue Cube
    elseif strcmp(cube_list(itr),'Blue')
        v1=[robot_b3d(1)-robot_r3d(1),robot_b3d(2)-robot_r3d(2)];
        v2=[cube_b(1)-robot_b3d(1),cube_b(2)-robot_b3d(2)];
        dot_product=dot(v1,v2);
        v1_norm=norm(v1);
        v2_norm=norm(v2);
        theta_b=acosd(dot_product/(v1_norm*v2_norm));
        cross_p=cross([v1,0], [v2,0]);
        if cross_p(3) > 0
            theta_b = -theta_b;
        end

        v1=[cube_b(1)-robot_b3d(1),cube_b(2)-robot_b3d(2)];
        v2=[circle_b(1)-cube_b(1),circle_b(2)-cube_b(2)];
        dot_product=dot(v1,v2);
        v1_norm=norm(v1);
        v2_norm=norm(v2);
        theta_cb=acosd(dot_product/(v1_norm*v2_norm));
        cross_p=cross([v1,0], [v2,0]);
        if cross_p(3) > 0
            theta_cb = -theta_cb;
        end
        d_tocube=sqrt((cube_b(1)-robot_b3d(1))^2+(cube_b(2)-robot_b3d(2))^2);
        d_tocircle=sqrt((cube_b(1)-circle_b(1))^2+(cube_b(2)-circle_b(2))^2);
        fprintf('turn(%f) ; go(%f) ; grab() ; turn(%f) ; go(%f) ; let_go() ; go(-10)\n',theta_b,d_tocube,theta_cb,d_tocircle-12)
        
%         d_10=(circle_b(1)-cube_b(1))^2+(circle_b(2)-cube_b(2))^2;
%         t_b=sqrt((10+13.5)^2/d_10);
%         t_r=sqrt((10+11.4+6)^2/d_10);
%         robot_b3d=[circle_b(1), circle_b(2)] + t_b * [cube_b(1)-circle_b(1), cube_b(2)-circle_b(2)];
%         robot_r3d=[circle_b(1), circle_b(2)] + t_r * [cube_b(1)-circle_b(1), cube_b(2)-circle_b(2)];
%         newp_b=M*[robot_b3d';7.5;1];
%         robot_b=newp_b./newp_b(3);
%         newp_r=M*[robot_r3d';7.5;1];
%         robot_r=newp_r./newp_r(3);

    end

end


%% Functions

function M = calibrate(points2d,points3d)

j = 1;

for i=1:length(points3d)
    A(j,:) = [points3d(1,i) points3d(2,i) points3d(3,i) 1 0 0 0 0 -points2d(1,i)*points3d(1,i) -points2d(1,i)*points3d(2,i) -points2d(1,i)*points3d(3,i) -points2d(1,i)];
    A(j+1,:) = [0 0 0 0 points3d(1,i) points3d(2,i) points3d(3,i) 1 -points2d(2,i)*points3d(1,i) -points2d(2,i)*points3d(2,i) -points2d(2,i)*points3d(3,i) -points2d(2,i)];

    j=j+2;
end

[~,~,V] = svd(A);

M = reshape(V(:,end),4,3);
end

function error=reprojection_error(M,p3d,p2d)
for i=1:size(p3d,2)
    pi=M*[p3d(:,i);1];
    pi=pi./pi(3);
    err(i)=norm(p2d(:,i)-pi(1:2));
end
error=sum(err)/length(err);
end

function [robot_b,robot_r] =  arrow(inputImage)
% Read the input image
inputImage = inputImage;

% Convert the input image to the HSV color space
hsvImage = rgb2hsv(inputImage);

% Extract the hue, saturation, and value channels
hueChannel = hsvImage(:, :, 1);
saturationChannel = hsvImage(:, :, 2);
valueChannel = hsvImage(:, :, 3);

% Define the range of values for a light white arrow in the HSV color space
hueThreshold = [0.05 0.2];
saturationThreshold = [0.0 0.2];
valueThreshold = [0.7 0.95];

% Create binary masks based on the threshold values
hueMask = (hueChannel >= hueThreshold(1)) & (hueChannel <= hueThreshold(2));
saturationMask = (saturationChannel >= saturationThreshold(1)) & (saturationChannel <= saturationThreshold(2));
valueMask = (valueChannel >= valueThreshold(1)) & (valueChannel <= valueThreshold(2));

% Combine the binary masks
binaryImage = hueMask & saturationMask & valueMask;

% Remove small objects and noise from the binary image
binaryImage = bwareaopen(binaryImage, 150);

% Perform morphological operations to enhance the arrow shape
se = strel('disk', 2);
morphedImage = imclose(binaryImage, se);

% Find arrow contours
contours = bwboundaries(morphedImage, 8, "noholes");

% Iterate through the detected contours and filter arrow shapes
for i = 1:length(contours)
    boundary = contours{i};
    perimeter = sum(sqrt(sum(diff(boundary).^2, 2)));
    area = polyarea(boundary(:, 2), boundary(:, 1));

    area_list(i) = area;
    if area < 800 || area > 3000
        continue;
    end

    % Calculate region properties including centroid
    regionProps = regionprops(morphedImage, 'Centroid', 'Area');

    % Get the centroid of the current region
    centroid = regionProps(i).Centroid;
    centroid_white_LED(i,:) = centroid;
    % Display the centroid on the original image with a marker
    inputImage = insertMarker(inputImage, centroid, 'o', 'Color', 'red', 'Size', 10);
end

% Display the original image with the centroids highlighted
imshow(inputImage);

% Convert image to HSV color space
hsvImage = rgb2hsv(inputImage);

% Define lower and upper range for blue color in HSV
lowerBlue = [0.40, 0.68, 0.3]; % Adjust these values as per your requirement
upperBlue = [0.55, 0.8, 0.45];

% Threshold the image to create a binary mask
blueMask = (hsvImage(:,:,1) >= lowerBlue(1) & hsvImage(:,:,1) <= upperBlue(1)) ...
    & (hsvImage(:,:,2) >= lowerBlue(2) & hsvImage(:,:,2) <= upperBlue(2)) ...
    & (hsvImage(:,:,3) >= lowerBlue(3) & hsvImage(:,:,3) <= upperBlue(3));

% Apply morphological operations (optional)
% se = strel('disk', 5);
% blueMask = imopen(blueMask, se);

% Find connected components in the mask
blueMask = bwlabel(blueMask);

blueMask= bwareaopen(blueMask, 5);

stats = regionprops(blueMask, 'BoundingBox', 'Centroid');

% Perform further analysis or actions based on the detected blue regions
imshow(inputImage);
hold on;
for i = 1:numel(stats)
    % Get the bounding box or centroid of the blue region
    bbox = stats(i).BoundingBox;
    centroid = stats(i).Centroid;
    centroid_blue_points(i,:) = centroid;

    % Draw the bounding box or centroid on the image (optional)
    rectangle('Position', bbox, 'EdgeColor', 'g', 'LineWidth', 2);
    plot(centroid(1), centroid(2), 'r+', 'MarkerSize', 10, 'LineWidth', 2);
end
hold off;

% Example points
vector1 = centroid_white_LED;
vector2 = centroid_blue_points;

% Calculate pairwise distances
distances = pdist2(vector1, vector2);

% Find the indices of the minimum distances
[minDist, minIndex] = min(distances(:));
[minRow, minCol] = ind2sub(size(distances), minIndex);

% Get the points with the minimum distance
point1 = vector1(minRow, :);
point2 = vector2(minCol, :);

% Display the results
disp("Points with the minimum distance:");
disp("Point 1: " + mat2str(point1));
disp("Point 2: " + mat2str(point2));
disp("Minimum distance: " + num2str(minDist));

% Calculate the middle point
middlePoint = (point1 + point2) / 2;

% Display the middle point
disp("Middle point: " + mat2str(middlePoint));


hold on
plot(point1(1), point1(2), 'b+', 'MarkerSize', 10, 'LineWidth', 2);
plot(point2(1), point2(2), 'r+', 'MarkerSize', 10, 'LineWidth', 2);
plot(middlePoint(1), middlePoint(2), 'g+', 'MarkerSize', 10, 'LineWidth', 2);
plot((point1(1)+middlePoint(1))/2,(point1(2)+middlePoint(2))/2,'c+','MarkerSize', 10, 'LineWidth', 2)

% robot_b=[(point1(1)+middlePoint(1))/2 (point1(2)+middlePoint(2))/2];
robot_b=[middlePoint(1) middlePoint(2)];
robot_r=[point1(1) point1(2)];

bluePoint = point1;
redPoint = point2;
COR = middlePoint;
end
