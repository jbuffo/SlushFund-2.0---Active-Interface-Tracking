function [T_source,S_source]=one_D_advective_flux2(phi_i,S_i,T_i,...
    ice_thickness,grid_size,dt,S_bottom,g,rho_sw)

%% Using formulation of Griewank and Notz (2013)
%% Constants
%g=10;                                % Gravity
%rho_sw=1028;                         % Density of sea water
%rho_sw=1012;
beta=0.0005836;                      % Coefficient for density dependence on salinity
h_i=[0:grid_size:ice_thickness];     % see GN (2013)
kappa=1.37*10^-7;                    % Thermal Diffusivity
mu=1.88*10^-3;                       % Viscosity
Ra_c=0.0101;                           % Critical Rayleigh #
%Ra_c=1.01;
S_sw=S_bottom;                             % Seawater salinity
alpha=1.56*10^-1;                    % Linear coeff for Rayleigh number driven advection
%alpha=1.56*10^1;
c_br=3985;                           % Specific heat of sea water

% Locate ice-ocean interface
depth=0;
for i=1:(ice_thickness/grid_size)+1;
    if phi_i(i)<1
        depth=i;
    else
        break;
    end
end

%% Conditions for advection
%% Local Rayleigh number calculation
pi_i=(10^-17)*((10^3)*phi_i).^3.1;          % permeability from GN (2013)
%pi_i=(5*10^-11)*((phi_i.^3)/(1-phi_i).^2);
% pi_i=0*[1:(ice_thickness/grid_size)+1];
% for i=1:(ice_thickness/grid_size)+1;
%     if phi_i(i)==1;
%         pi_i(i)=1;
%     else
%         pi_i(i)=heaviside(0.849734-phi_i(i)).*(0.5+(1/pi)*atan(10*(0.5-phi_i(i)))).*...
%             (5*10^-9).*((phi_i(i).^3)./(1-phi_i(i)).^2)+heaviside(phi_i(i)-0.849734).*...
%             (10^-17).*((10^3)*phi_i(i)).^3.1;
%     end
% end
% for i=1:length(phi_i);
%     if phi_i(i)<=0.05;
%         pi_i(i)=0;
%     else
%         pi_i(i)=pi_i(i);
%     end
% end
Ra_i=0*[1:(ice_thickness/grid_size)+1];
for i=1:(ice_thickness/grid_size)+1;
    Ra_i(i)=(g*rho_sw*beta*(S_i(i)-S_sw)*min(pi_i(i:end))*(grid_size*depth-h_i(i)))/(kappa*mu);    % Local Rayleigh #
end

%% Calculating fluxes
% Flux into brine channel
F_c=0*[1:(ice_thickness/grid_size)+1];
for i=1:(ice_thickness/grid_size)+1;
    if i>depth;
        F_c(i)=0;
    elseif Ra_i(i)<Ra_c;
        F_c(i)=0;
    elseif phi_i(i)<0.05;
        F_c(1:i)=0;
    else
        F_c(i)=alpha*(Ra_i(i)-Ra_c)*(grid_size^3)*dt;
    end
end
% Vertical Flux into and out of cell
F_in=0*[1:(ice_thickness/grid_size)+1];
F_out=0*[1:(ice_thickness/grid_size)+1];
for i=1:(ice_thickness/grid_size)+1;
    if i==1;
        F_out(i)=0;
        if F_c(i)==0
            F_in(i)=0;
        else
            F_in(i)=F_c(i);
        end
    elseif phi_i(i)>=1
        F_in(i)=0;
        F_out(i)=0;
    else
        F_out(i)=F_in(i-1);
        F_in(i)=F_c(i)+F_out(i);
    end
end

%% Calculate sources or sinks for Temperature and Salinity
T_source=0*[1:(ice_thickness/grid_size)+1];
S_source=0*[1:(ice_thickness/grid_size)+1];
for i=1:(ice_thickness/grid_size)+1;
    if i==1
        S_source(i)=(S_i(i+1)*F_in(i))-(S_i(i)*F_c(i));
        T_source(i)=c_br*(rho_sw+rho_sw*beta*(S_i(i)-S_sw))...
            *(grid_size^3)*(((F_in(i)*T_i(i+1)+(((rho_sw+rho_sw*beta*(S_i(i)-S_sw))...
            *(grid_size^3))-F_in(i))*T_i(i))/...
            ((rho_sw+rho_sw*beta*(S_i(i)-S_sw))...
            *(grid_size^3)))-T_i(i));
    elseif i<(ice_thickness/grid_size)+1
        S_source(i)=(S_i(i+1)*F_in(i))-(S_i(i)*F_c(i))-...
            (S_i(i)*F_out(i));
        T_source(i)=c_br*(rho_sw+rho_sw*beta*(S_i(i)-S_sw))...
            *(grid_size^3)*(((F_in(i)*T_i(i+1)+(((rho_sw+rho_sw*beta*(S_i(i)-S_sw))...
            *(grid_size^3))-F_in(i))*T_i(i))/...
            ((rho_sw+rho_sw*beta*(S_i(i)-S_sw))...
            *(grid_size^3)))-T_i(i));
    else
        S_source(i)=0;
        T_source(i)=0;
    end
end

T_source=T_source';
S_source=S_source';

