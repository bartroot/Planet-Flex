function [M] = MaternCovariance_2(x, sigma, nu, rho, tf_dist)
%MaternCovariance(x1, x2, sigma, nu, alpha): Function that calculates
%covariance matrix with covariances between points in array x based on the 
%distance between the points, according to the Mátern class of covariance 
%function 
% 
% Output:   M       Covariance Matrix
% 
% Input:    x       [m x 3]-matrix with points per row and coordinates over 
%                   columns
%           sigma   standard deviation
%           nu      smoothness parameter
%           alpha   scale paramter
% 
% Reference article: 
%   Matérn, B. (1960). Spatial variation, meddelanden fran statens 
%       skogsforskningsinstitut. Lecture Notes in Statistics, 36:21.
m=length(x);
x(:,1) = x(:,1)-1;
M=eye(m)*sigma^2;
distance = eye(m);

if tf_dist == 1
for i= 1:m
    j=i+1;
    while j<=m
        
        if x(i,1)>= 0 && x(j,1)>=0
         dist1 = (x(i,1)-x(j,1)).^2;
        elseif x(i,1)<= 0 && x(j,1)>=0
         dist1 = (x(i,1)-x(j,1)).^2;
         if sqrt(dist1) > 180
             dist1 = (360 - sqrt(dist1)).^2;
         end
        elseif x(i,1)>= 0 && x(j,1)<=0
         dist1 = (x(i,1)-x(j,1)).^2;
         if sqrt(dist1) > 180
             dist1 = (360 - sqrt(dist1)).^2;
         end
        elseif x(i,1)<= 0 && x(j,1)<=0
         dist1 = (x(i,1)-x(j,1)).^2;
        end
        
        distance(i,j)  = sind(real(acosd(sind(x(i,2))*sind(x(j,2)) + cosd(x(i,2))*cosd(x(j,2))*cosd(sqrt(dist1))))/2);
        distance(j,i) = distance(i,j);
        j=j+1;
    end
end
    save( 'Functions_Heterogeneity/distance_matrix.mat', 'distance' );

else
    load('Functions_Heterogeneity/distance_matrix.mat');
end

breaks = 2;
for i = 1:breaks
    j = (i-1)*(m/breaks)+1;
    k = (i)*(m/breaks);
    test = sigma.^2.*2.^(1-nu)./gamma(nu).*(4*sqrt(nu)/rho.*distance(j:k,j:k)).^nu.*besselk(nu,4*sqrt(nu)/rho.*distance(j:k,j:k));
    M(j:k,j:k) = test-(test(1,1)*eye(m/100))+(sigma^2*eye(m/100));
end
        %M = sigma.^2.*2.^(1-nu)./gamma(nu).*(4*sqrt(nu)/rho.*sind(distance./2)).^nu.*besselk(nu,4*sqrt(nu)/rho.*sind(distance./2));
end

