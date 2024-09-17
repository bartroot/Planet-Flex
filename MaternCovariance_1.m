function [M] = MaternCovariance_1(x, sigma, nu, rho)
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
for i= 1:m
    j=i+1;
    while j<=m
        
        if x(i,1)>= 0 && x(j,1)>=0
         dist1 = (x(i,1)-x(j,1)).^2;
        elseif x(i,1)<= 0 && x(j,1)>=0
         dist1 = (x(i,1)-x(j,1)).^2;
         if sqrt(dist1) >= 180
             dist1 = (360 - sqrt(dist1)).^2;
         end
        elseif x(i,1)>= 0 && x(j,1)<=0
         dist1 = (x(i,1)-x(j,1)).^2;
         if sqrt(dist1) >= 180
             dist1 = (360 - sqrt(dist1)).^2;
         end
        elseif x(i,1)<= 0 && x(j,1)<=0
         dist1 = (x(i,1)-x(j,1)).^2;
        end
        
        distance  = acosd(sind(x(i,2))*sind(x(j,2)) + cosd(x(i,2))*cosd(x(j,2))*cosd(sqrt(dist1)));
        
        M(i,j) = sigma.^2.*2.^(1-nu)./gamma(nu).*(4*sqrt(nu)/rho.*sind(distance/2)).^nu.*besselk(nu,4*sqrt(nu)/rho.*sind(distance/2));
        M(j,i)=M(i,j);
        j=j+1;
    end
end
end

