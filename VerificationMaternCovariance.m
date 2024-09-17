% Thesis Stan Spee 2021
% Verification of Matérn covariance function
% Ref Article: 
% François Bachoc. Parametric estimation of covariance function in Gaussian-process based Kriging
% models. Application to uncertainty quantification for computer experiments. Statistics [math.ST].
% Université Paris-Diderot - Paris VII, 2013. English. fftel-00881002f
% And:
% Thesis of Robin Thor
close all;
fig=1;
%% Dependency on nu
%values by thor:
x=[0:0.001:deg2rad(50)]';
x=[x,zeros(length(x),2)];
nu=[.1,.3,.5,.7,.9,1.1,1.3,1.5,1.7];
l=2;
sigma = sqrt(10);
for i=1:length(nu)
    CovarianceMatrix_vv=MaternCovariance(x, sigma, nu(i), 2*sqrt(nu(i))/l);
    figure(fig);
    plot(rad2deg(x(:,1)),CovarianceMatrix_vv(1,:),'LineWidth',1.5);hold on
    xlabel("\psi (deg)");ylabel("Covariance (km^2)");
    title(["\it{l}="+num2str(l)+"; \sigma="+num2str(sigma)])

    figure(fig+1)
    rng(1)
    y=mvnrnd(zeros(length(x),1), CovarianceMatrix_vv); 
    plot(rad2deg(x(:,1)),rad2deg(y),'LineWidth',1.5);hold on
    xlabel("\psi (deg)");ylabel("y");
    title(["\it{l}="+num2str(l)+"; \sigma="+num2str(sigma)])
    ind(i)=["\nu = "+ num2str(nu(i))];
end
figure(fig);legend(ind);
figure(fig+1);legend(ind);
fig=fig+2;
%values by Bachoc:
x=[0:0.01:1.9]';
x=[x,zeros(length(x),2)];
nu=[0.01 0.5 1.5 2.5 100];
l=1;
sigma=300;
for i=1:length(nu)
    CovarianceMatrix_vv=MaternCovariance(x, sigma, nu(i), 2*sqrt(nu(i))/l);
    figure(fig)
    plot(x(:,1),CovarianceMatrix_vv(1,:),'LineWidth',1.5);hold on
    xlabel("Distance (m)");ylabel("Covariance (m^2)");
    title(["\it{l}="+num2str(l)+"; \sigma="+num2str(sigma)])
    
    figure(fig+1)
    rng(1)
    y=mvnrnd(zeros(length(x),1), CovarianceMatrix_vv); 
    plot(x(:,1),y,'LineWidth',1.5);hold on
    xlabel("Distance (m)");ylabel("y");
    title(["\it{l}="+num2str(l)+"; \sigma="+num2str(sigma)])
    ind(i)=["\nu = "+ num2str(nu(i))];
end
figure(fig);legend(ind)
figure(fig+1);legend(ind)
fig=fig+2;

%% Dependency on l
x=[0:0.001:deg2rad(25)]';
x=[x,zeros(length(x),2)];
l=[.01,.032,.1,.32,1,3.2,10];
nu=1;
sigma = sqrt(10);
for i=1:length(l)
    CovarianceMatrix_vv=MaternCovariance(x, sigma, nu, 2*sqrt(nu)/l(i));
    figure(fig);
    plot(rad2deg(x(:,1)),CovarianceMatrix_vv(1,:),'LineWidth',1.5);hold on
    xlabel("\psi (deg)");ylabel("Covariance (km^2)");
    title(["\nu="+num2str(nu)+"; \sigma="+num2str(sigma)])

    figure(fig+1)
    rng(2)
    y=mvnrnd(zeros(length(x),1), CovarianceMatrix_vv); 
    plot(rad2deg(x(:,1)),rad2deg(y),'LineWidth',1.5);hold on
    xlabel("\psi (deg)");ylabel("y");
    title(["\nu="+num2str(nu)+"; \sigma="+num2str(sigma)])
    ind(i)=["\it{l} = "+ num2str(l(i))];
end
figure(fig);legend(ind);
figure(fig+1);legend(ind);
fig=fig+2;
%values by Bachoc:
x=[0:0.01:2]';
x=[x,zeros(length(x),2)];
l=[0.5 1.5 2.5 4 5 ];
nu=1;
sigma=1;
for i=1:length(l)
    CovarianceMatrix_vv=MaternCovariance(x, sigma, nu, 2*sqrt(nu)/l(i));
    figure(fig)
    plot(x(:,1),CovarianceMatrix_vv(1,:),'LineWidth',1.5);hold on
    xlabel("Distance (m)");ylabel("Covariance (m^2)");
    title(["\nu="+num2str(nu)+"; \sigma="+num2str(sigma)])
    
    figure(fig+1)
    rng(2)
    y=mvnrnd(zeros(length(x),1), CovarianceMatrix_vv); 
    plot(x(:,1),y,'LineWidth',1.5);hold on
    xlabel("Distance (m)");ylabel("y");
    title(["\nu="+num2str(nu)+"; \sigma="+num2str(sigma)])
    ind(i)=["\it{l} = "+ num2str(l(i))];
end
figure(fig);legend(ind)
figure(fig+1);legend(ind)
fig=fig+2;

%% Dependency on \sigma
x=[0:0.0001:deg2rad(25)]';
x=[x,zeros(length(x),2)];
l=1;
nu=1;
sigma = sqrt([1 5 10 20]);
for i=1:length(sigma)
    CovarianceMatrix_vv=MaternCovariance(x, sigma(i), nu, 2*sqrt(nu)/l);
    figure(fig);
    plot(rad2deg(x(:,1)),CovarianceMatrix_vv(1,:),'LineWidth',1.5);hold on
    xlabel("\psi (deg)");ylabel("Covariance (km^2)");
    title(["\nu="+num2str(nu)+"; \it{l}="+num2str(l)])

    figure(fig+1)
    rng(2)
    y=mvnrnd(zeros(length(x),1), CovarianceMatrix_vv); 
    plot(rad2deg(x(:,1)),rad2deg(y),'LineWidth',1.5);hold on
    xlabel("\psi (deg)");ylabel("y");
    title(["\nu="+num2str(nu)+"; \it{l}="+num2str(l)])
    ind(i)=["\sigma = "+ num2str(sigma(i))];
end
figure(fig);legend(ind);
figure(fig+1);legend(ind);
fig=fig+2;
%values by Bachoc:
x=[0:0.01:2]';
x=[x,zeros(length(x),2)];
l= 1;
nu=1;
sigma= [0.5 1 1.5 2];
for i=1:length(sigma)
    CovarianceMatrix_vv=MaternCovariance(x, sigma(i), nu, 2*sqrt(nu)/l);
    figure(fig)
    plot(x(:,1),CovarianceMatrix_vv(1,:),'LineWidth',1.5);hold on
    xlabel("Distance (m)");ylabel("Covariance (m^2)");
    title(["\nu="+num2str(nu)+"; \it{l}="+num2str(l)])
    
    figure(fig+1)
    rng(2)
    y=mvnrnd(zeros(length(x),1), CovarianceMatrix_vv); 
    plot(x(:,1),y,'LineWidth',1.5);hold on
    xlabel("Distance (m)");ylabel("y");
    title(["\nu="+num2str(nu)+"; \it{l}="+num2str(l)])
    ind(i)=["\sigma = "+ num2str(sigma(i))];
end
figure(fig);legend(ind)
figure(fig+1);legend(ind)
fig=fig+2;





%%
CenterOfMass=asteroid_ht.Polyhedron.CentersOfMass;
x=CenterOfMass./max(CenterOfMass);
distance=sqrt(sum((x(:,:)-x(1,:)).^2,2));%distance of element 1 with other elements(for figure of covariance)
nu=[0.5 1.5 2.5 100];
l=1;
sigma=1;
for i=1:length(nu)
    CovarianceMatrix_vv=MaternCovariance(x, sigma, nu(i), 2*sqrt(nu(i))/l);
    figure(fig)
    scatter(distance,CovarianceMatrix_vv(1,:),'LineWidth',1.5);hold on
    xlabel("Distance (m)");ylabel("Covariance (m^2)");
    title(["\it{l}="+num2str(l)+"; \sigma="+num2str(sigma)])
    
    figure(fig+1)
    rng(1)
    y=mvnrnd(zeros(length(x),1), CovarianceMatrix_vv); 
    scatter(distance,y,'LineWidth',1.5);hold on
    xlabel("Distance (m)");ylabel("y");
    title(["\it{l}="+num2str(l)+"; \sigma="+num2str(sigma)])
    ind(i)=["\nu = "+ num2str(nu(i))];
end
figure(fig);legend(ind)
figure(fig+1);legend(ind)


