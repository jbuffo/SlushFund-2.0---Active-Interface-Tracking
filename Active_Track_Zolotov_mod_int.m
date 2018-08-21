function history=Active_Track_Zolotov_mod_int(TBOTTOM, SBOTTOM,...
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

%T_grad=0.080766;
%T_grad=(Tbottom-Ttop)/Interface_depth;

%% BC=1 for Dirilecht Boundary Condition, BC=2 Neumann BC
%BC=1;
BC=2;

Depth=[0:-1:-(H/dz)];
%Precipitates=zeros((H/dz)+1,50);
%Brine_comp=zeros((H/dz)+1,15);
%for j=1:(H/dz)+1
%    Brine_comp(j,:)=[0 1 1 0.469 0.0102 0.0103 0.0528 0.546 0.0282 0.0 0.0 0.0...
%     273.15 272.15 1.0];
%end

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
%Thold=[];
Liquid_Fraction=[];
%phihold=[];
Salinity=[];
%Shold=[];
Darcy_Velocity=[];
%whold=[];
%Temp_mov=[];


%% Looping Over Time
for n=1:tf/dt
    if n==1;
        %[NewTemp,NewPhi,NewS,Neww]=...
        %    ESM_project_T_phi_S_advection_bars2(k_i,k_br,...
        %    k_s,dt,dz,H,c_i,c_br,L,Tol,phi_initial,T_initial,...
        %    S_initial,S_bottom,rho_i,rho_br,Tm,m,w_initial,n,Ttop,Tbottom);
        [NewTemp,NewPhi,NewS,Neww]=one_D_adv_ARD_FREZ_track(...
        k_i,k_br,k_s,dt,dz,H,c_i,c_br,L,STol,TTol,PhiTol,phi_initial,...
        T_initial,S_initial,S_bottom,rho_i,...
        rho_br,Tm,m,w_initial,n,Ttop,Tbottom,T_grad,BC,g,rho_br);
        %NewTemp(length(NewTemp))=Tbottom;
      
        %for j=1:(H/dz)+1
        %    Brine_comp(j,4:9)=(NewS(j)/S_initial(j))*Brine_comp(j,4:9);
        %    Brine_comp(j,13)=T_initial(j);
        %    Brine_comp(j,14)=NewTemp(j);
        %    Brine_comp(j,15)=abs(T_initial(j)-NewTemp(j));
        %    S_in=sum(Brine_comp(j,4:9));
        %    [conditions,solution,solid]=FREZCHEM(Brine_comp(j,:));
        %    Brine_comp(j,4:9)=solution(:,5);
        %    S_out=sum(Brine_comp(j,4:9));
        %    NewS(j)=(S_out/S_in)*NewS(j);
        %    Precipitates(j,:)=Precipitates(j,:)+(NewPhi(j)*solid(:,1))';
        %end
            
        %NewS(1)=0;
        %NewS(length(NewS))=S_bottom;
        Temperature=NewTemp;
        %Thold=[Thold NewTemp];
        Liquid_Fraction=NewPhi;
        %phihold=[phihold NewPhi];
        Salinity=NewS;
        %Shold=[Shold NewS];
        Darcy_Velocity=Neww;
        %whold=[whold Neww];
    else
        [NewTemp,NewPhi,NewS,Neww]=one_D_adv_ARD_FREZ_track(...
        k_i,k_br,k_s,dt,dz,H,c_i,c_br,L,STol,TTol,PhiTol,Liquid_Fraction',...
        Temperature',Salinity',S_bottom,rho_i,...
        rho_br,Tm,m,Darcy_Velocity',n,Ttop,Tbottom,T_grad,BC,g,rho_br);
        %NewTemp(1)=Ttop;
        %NewTemp(length(NewTemp))=Tbottom;
        
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
        
        %for j=1:(H/dz)+1
        %    Brine_comp(j,4:9)=(NewS(j)/S_initial(j))*Brine_comp(j,4:9);
        %    Brine_comp(j,13)=T_initial(j);
        %    Brine_comp(j,14)=NewTemp(j);
        %    Brine_comp(j,15)=abs(T_initial(j)-NewTemp(j));
        %    S_in=sum(Brine_comp(j,4:9));
        %    [conditions,solution,solid]=FREZCHEM(Brine_comp(j,:));
        %    Brine_comp(j,4:9)=solution(:,5);
        %    S_out=sum(Brine_comp(j,4:9));
        %    NewS(j)=(S_out/S_in)*NewS(j);
        %    Precipitates(j,:)=Precipitates(j,:)+(NewPhi(j)*solid(:,1))';
        %end
        
        %NewS(1)=0;
        %NewS(length(NewS))=S_bottom;
        
        %for i=1:(H/dz)+1;
        %    if NewPhi(i)>0.9;
        %        NewTemp(i)=Tm-(1.853*(34/28))-0.05;
        %        NewS(i)=34;
        %    else
        %        NewTemp(i)=NewTemp(i);
        %    end
        %end
        %Temperature=NewTemp;
        %Thold=[Thold NewTemp];
        %Liquid_Fraction=NewPhi;
        %phihold=[phihold NewPhi];
        %Salinity=NewS;
        %Shold=[Shold NewS];
        %Darcy_Velocity=Neww;
        %whold=[whold Neww];


%         if n==500 || n==1000 || n==1500 || n==2000 || n==4000 || n==6000 || n==8000 || n==10000 || n==12000 ...
%             || n==14000 || n==16000 || n==18000;
%             disp(n)
%             datestr(now)
%         else
%         end
        %% Movie creation (comment out unless making movies)
        %set(gcf,'visible','off')
        %fig=figure;
        %    subplot(1,5,1);
        %    image(NewTemp,'CDataMapping','scaled');
        %    title('Temperature');
        %    colorbar;
        %    subplot(1,5,2);
        %    image(NewPhi,'CDataMapping','scaled');
        %    title('Liquid Fraction');
        %    colorbar;
        %    subplot(1,5,3);
        %    image(NewS,'CDataMapping','scaled');
        %    title('Salinity');
        %    colorbar;
        %    subplot(1,5,4);
        %    image(Neww,'CDataMapping','scaled');
        %    title('Brine Velocity');
        %    colorbar;
        %    subplot(1,5,5);
        %    image(Salt_in_Ice(:,n-1)+IceSalt,'CDataMapping','scaled');
        %    title('Salt in Ice');
        %    colorbar;
        %temp=getframe(fig);
        %Temp_mov=[Temp_mov temp];
        %% End of movie making code
        
    end
end

Tplot=[Tplot Temperature];
Splot=[Splot Salinity];
phiplot=[phiplot Liquid_Fraction];
wplot=[wplot Darcy_Velocity];
Stotplot=[Stotplot Liquid_Fraction.*Salinity];


%% Output Images
%figure
%subplot(2,4,1)
%image(Temperature,'CDataMapping','scaled')
%title('Temperature')
%xlabel(['Time (' num2str(dt) 'sec)'])
%ylabel('Depth (cm)')
%colorbar
%subplot(2,4,2)
%image(Liquid_Fraction,'CDataMapping','scaled')
%title('Liquid Fraction')
%xlabel(['Time (' num2str(dt) 'sec)'])
%ylabel('Depth (cm)')
%colorbar
%subplot(2,4,3)
%image(Salinity,'CDataMapping','scaled')
%title('Salinity')
%xlabel(['Time (' num2str(dt) 'sec)'])
%ylabel('Depth (cm)')
%colorbar
%subplot(2,4,4)
%image(Darcy_Velocity(1:end,10:end),'CDataMapping','scaled')
%title('Darcy Velocity')
%xlabel(['Time (' num2str(dt) 'sec)'])
%ylabel('Depth (cm)')
%colorbar
%subplot(2,4,5)
%image(Liquid_Fraction.*Salinity,'CDataMapping','scaled')
%title('Total Salt (ppt)')
%xlabel(['Time (' num2str(dt) 'sec)'])
%ylabel('Depth (cm)')
%colorbar
%subplot(2,4,6)
%image(Platelet_frac,'CDataMapping','scaled')
%title('Platelet Ice Fraction')
%xlabel(['Time (' num2str(dt) 'sec)'])
%ylabel('Depth (cm)')
%colorbar
%subplot(2,4,7)
%plot(Platelet_frac(:,end),Depth,0*Depth+0.1,Depth)
%title('Platelet Fraction vs. Depth')
%xlabel('Platelet Fraction')
%ylabel('Depth (cm)')
%subplot(2,4,8)
%text(0.0,0.0,['Supercooled ' num2str(supercooling) 'K'])
%text(0.0,0.1,['Ttop= ' num2str(Ttop) ' Kelvin'])
%text(0.0,0.2,['Tbottom= ' num2str(Tbottom) ' Kelvin'])
%text(0.0,0.3,['tf= ' num2str(t_end/86400) ' days'])
%text(0.0,0.4,['K= ' num2str(K)])
%text(0.0,0.5,['dt= ' num2str(dt) ' sec'])
%axis off

%implay(Temp_mov)
%movie2avi(Temp_mov,'ice_shelf.avi')

save(strcat('HECCtest_Europa_',num2str(interface_in),'m_salinity_',num2str(S_bottom)...
    ,'_Tgrad_',num2str(T_grad),'_time_',num2str(t_start),'_to_',num2str(t_end),...
    '.mat'));

%% Initial Conditions for next time loop
% IC=[Temperature Salinity Liquid_Fraction Darcy_Velocity];
% end
% 
% figure
% subplot(2,4,1)
% plot(Tplot,Depth)
% title('Temperature vs. Depth over time')
% xlabel('Temperature (K)')
% ylabel('Depth (cm)')
% subplot(2,4,2)
% plot(phiplot,Depth)
% title('Liquid Fraction vs. Depth over time')
% xlabel('Liquid Fraction')
% ylabel('Depth (cm)')
% subplot(2,4,3)
% plot(Splot,Depth)
% title('Salinity vs. Depth over time')
% xlabel('Salinity (psu)')
% ylabel('Depth (cm)')
% subplot(2,4,4)
% plot(wplot,Depth)
% title('Darcy Velocity vs. Depth over time')
% xlabel('Darcy Velocity (m/s)')
% ylabel('Depth (cm)')
% subplot(2,4,5)
% plot(phiplot.*Splot,Depth)
% title('Total Salt vs. Depth over time')
% xlabel('Total Salt (ppt)')
% ylabel('Depth (cm)')
% subplot(2,4,8)
% %text(0.0,0.1,['Ttop= ' num2str(Ttop) ' Kelvin'])
% %text(0.0,0.2,['Tbottom= ' num2str(Tbottom) ' Kelvin'])
% text(0.0,0.3,['tf= ' num2str(t_end/86400) ' days'])
% %text(0.0,0.4,['K= ' num2str(K)])
% text(0.0,0.5,['dt= ' num2str(dt) ' sec'])
% axis off

[Temperature Salinity Liquid_Fraction]
[history]