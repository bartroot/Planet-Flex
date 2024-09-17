function [M] = MaternCovariance(x, sigma, nu, alpha)
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
M=eye(m)*sigma^2;
for i= 1:m
    j=i+1;
    while j<=m
        distance=sqrt(sum((x(i,:)-x(j,:)).^2,2));
        M(i,j) = sigma.^2.*2.^(1-nu)./gamma(nu).*(alpha.*distance).^nu.*besselk(nu,alpha.*distance);
        M(j,i)=M(i,j);
        j=j+1;
    end
end
end

