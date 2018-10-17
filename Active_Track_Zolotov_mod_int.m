function [Temperature, Salinity, Liquid_Fraction, history]=...
    Active_Track_Zolotov_mod_int(TBOTTOM, SBOTTOM,...
    INTERFACE_DEPTH, T_END, RHO_BR, DT)

%% Vectors for each time chunk of iteration for plotting
Tplot=[];
Splot=[];
phiplot=[];
wplot=[];
Stotplot=[];

%% Time loops
%time=[T_START:T_END:T_END];
%for k=1:length(time)-1;

%% Input Parameters 
Ttop=100;                %% Ice shelf top temp
%Tbottom=273.0;
Tbottom=TBOTTOM;
%interface_in=10;         %% Initial interface depth (m)
interface_in=INTERFACE_DEPTH;
%Tbottom=271.2;             %% '       ' bottom temp (-1.95C)
Tm=273.15;               %% Melt temperature of pure ice
k_i=2;                   %% Diffusion Coefficient for Heat (ice)
k_br=.6;                 %% '                            ' (brine)
k_s=2*10^-9;                   %% Diffusion Coefficient for Salt
%dt=5000;                    %% Time step
dt=DT;
dz=.01;                   %% Spatial step
H=1.00;                   %% Thickness of ice shelf
c_br=3985;               %% Specific heat of seawater
c_i=2000;                %% '              ' ice
L=334774;                %% Latent heat of fusion ice<->water
%S_bottom=12.3;
S_bottom=SBOTTOM;
%S_bottom=34;             %% '         ' bottom of shelf
rho_i=917;               %% Density of Ice
%rho_br=1012;
rho_br=RHO_BR;
%rho_br=1025;             %% Density of Brine
% t_start=time(k);
% t_end=time(k+1);
% tf=t_end-t_start;                %% Final time
t_start=0;
t_end=T_END;
tf=T_END;
STol=.05;                  %% Error Tolerance (includes T&Phi&S)
PhiTol=0.005;
TTol=0.05;
m=2;                     %% Cementation exponent
%g=9.8;                    %% Earth Gravity
g=1.32;
mu=1.88*10^-3;           %% Viscosity
phi_c=0.05;              %% Critical Porosity


%% BC=1 for Dirilecht Boundary Condition, BC=2 Neumann BC
%BC=1;
BC=2;

Depth=[0:-1:-(H/dz)];


%% Initial Condition Vectors for Phi and T and w and S

if t_start==0
    T_initial=0*[0:dz:H]+Tbottom;
    S_initial=0*[0:dz:H]+S_bottom;
    phi_initial=0*[0:dz:H]+1;
    w_initial=0*[0:dz:H];
    history=[];
    Interface_depth=interface_in;
    T_grad=(Tbottom-Ttop)/Interface_depth;
else
    T_initial=IC(:,1)';
    S_initial=IC(:,2)';
    phi_initial=IC(:,3)';
    w_initial=IC(:,4)';
end

%% Vectors for plotting
Temperature=[];

Liquid_Fraction=[];

Salinity=[];

Darcy_Velocity=[];



%% Looping Over Time
for n=1:tf/dt
    if n==1;

        [NewTemp,NewPhi,NewS,Neww]=one_D_adv_ARD_FREZ_track(...
        k_i,k_br,k_s,dt,dz,H,c_i,c_br,L,STol,TTol,PhiTol,phi_initial,...
        T_initial,S_initial,S_bottom,rho_i,...
        rho_br,Tm,m,w_initial,n,Ttop,Tbottom,T_grad,BC,g,rho_br);

            

        Temperature=NewTemp;

        Liquid_Fraction=NewPhi;

        Salinity=NewS;

        Darcy_Velocity=Neww;

    else
        [NewTemp,NewPhi,NewS,Neww]=one_D_adv_ARD_FREZ_track(...
        k_i,k_br,k_s,dt,dz,H,c_i,c_br,L,STol,TTol,PhiTol,Liquid_Fraction',...
        Temperature',Salinity',S_bottom,rho_i,...
        rho_br,Tm,m,Darcy_Velocity',n,Ttop,Tbottom,T_grad,BC,g,rho_br);
        
        crit_interface=0;
        for i=1:H/dz+1
            if NewPhi(i)<=0.05;
                crit_interface=crit_interface+1;
            else
            break
            end
        end
        
        if crit_interface==0;
            Temperature=NewTemp;
            Liquid_Fraction=NewPhi;
            Salinity=NewS;
            Darcy_Velocity=Neww;
        else
            Temperature=[NewTemp(crit_interface+1:end); Tbottom*(1+0*[1:crit_interface])'];
            Liquid_Fraction=[NewPhi(crit_interface+1:end); 1*(1+0*[1:crit_interface])'];
            Salinity=[NewS(crit_interface+1:end); S_bottom*(1+0*[1:crit_interface])'];
            Darcy_Velocity=[Neww(crit_interface+1:end); Neww(end)*(1+0*[1:crit_interface])'];
            
            history=[history; NewTemp(1:crit_interface) NewS(1:crit_interface) ...
                NewPhi(1:crit_interface) Neww(1:crit_interface)];
            
            Interface_depth=Interface_depth+dz*crit_interface;
            T_grad=(Tbottom-Ttop)/Interface_depth;
        end
        
    end
end

Tplot=[Tplot Temperature];
Splot=[Splot Salinity];
phiplot=[phiplot Liquid_Fraction];
wplot=[wplot Darcy_Velocity];
Stotplot=[Stotplot Liquid_Fraction.*Salinity];


save(strcat('HECCtest_Europa_',num2str(interface_in),'m_salinity_',num2str(S_bottom)...
    ,'_Tgrad_',num2str(T_grad),'_time_',num2str(t_start),'_to_',num2str(t_end),...
    '.mat'));


[Temperature Salinity Liquid_Fraction]
[history]