function [T_np1_k_j,phi_np1_k_j,S_np1_k_j,w_np1_k_j]=one_D_adv_ARD_FREZ_track(...
    k_i,k_br,k_s,dt,dz,H,c_i,c_br,L,STol,TTol,PhiTol,phi_initial,T_initial,...
    S_initial,S_bottom,rho_i,rho_br,Tm,m,w_initial,n,Ttop,Tbottom,T_grad,BC,g,rho_sw)

%% Initial Condition Vectors for Phi, T, S, w
phi_initial=phi_initial';
T_initial=T_initial';
phi_np1_km2_j=phi_initial;
T_np1_km1_j=T_initial;
S_initial=S_initial';
S_np1_km1_j=S_initial;
w_initial=w_initial';
w_np1_k_j_loop=w_initial;

S_evolve=[]';
T_evolve=[]';
phi_evolve=[]';
SErr=STol+1;
TErr=TTol+1;
PhiErr=PhiTol+1;
counter=0;
matrix_dimension=H/dz+1;


%% Itterattions over k
while SErr>STol || TErr>TTol || PhiErr>PhiTol;
    counter=counter+1;
    if counter>25
        disp(counter)
    else
    end
    
    [T_source,S_source]=one_D_advective_flux2(phi_np1_km2_j,S_np1_km1_j...
    ,T_np1_km1_j,H,dz,dt,S_bottom,g,rho_sw);
    
%   delta_T=1.853*(S_np1_km1_j/28);        %% Old Freezing point depression due to salt
    delta_T=Tm-(-(1.333489497*(10^-5)*S_np1_km1_j.^2)-0.01612951864*S_np1_km1_j+...
        273.055175687);                       %% Liquidus point using FREZCHEM Europa ocean
%   delta_T=Tm-(-(9.1969758*(10^-5)*S_np1_km1_j.^2)-0.03942059*S_np1_km1_j+...
%        272.63617665);                        %% Liquidus point using FREZCHEM seawater
    Hs=c_i*(Tm-delta_T);                     %% Enthalpy Calculation
    Hsmelt=c_i*Tm;                           %% Enthalpy to melt ice
    %% Calculating solid fraction at k-1 iteration for use in
    %% temperature, salinity, and brine velocity at k iteration
    
    %% Value of Phi(n+1,k-1,j)
    phi_np1_km1_j=phi_np1_km2_j;
    for i=1:matrix_dimension;
    %% Test to see if melting ice or freezing brine
    En=c_i*T_np1_km1_j(i)+L*phi_np1_km2_j(i);
        if En<Hs(i)
            phi_np1_km1_j(i)=0;
        elseif En>Hs(i)+L;
            phi_np1_km1_j(i)=1;
        else
            phi_np1_km1_j(i)=((c_i*T_np1_km1_j(i)+phi_np1_km2_j(i)*L)-Hs(i))/...
            L;
        end;
    end;

    %% Reassigning Phi(n+1,k-2) as Phi(n+1,k-1) for next
    %% iteration
    phi_np1_km2_j=phi_np1_km1_j;
    
    %% Solving for Darcy Velocity (eq. 9 Hunke et. al.) no gravity drainage

         w_np1_k_j=w_initial;

    w_np1_k_j=w_np1_k_j_loop;
    w_np1_k_j_loop=w_initial;
    
    %% Create K matrix and solve for T
    D=((phi_np1_km1_j*k_br+(1-phi_np1_km1_j)*k_i)*dt)./...
        ((dz^2)*(phi_np1_km1_j*(rho_br*c_br)+(1-phi_np1_km1_j)...
        *(rho_i*c_i)));
    D_w=rho_br*c_br*w_np1_k_j*dt./(dz*(phi_np1_km1_j*(rho_br*c_br)+(1-phi_np1_km1_j)...
        *(rho_i*c_i)));
    for i=1:matrix_dimension;
        if w_np1_k_j(i)==0;
            a(i)=1+2*D(i);
            b(i)=-D(i);
            c(i)=-D(i);
        elseif w_np1_k_j(i)>0;
            a(i)=1+2*D(i)+D_w(i);
            b(i)=-D(i);
            c(i)=-D(i)-D_w(i);
        else a(i)=1+2*D(i)-D_w(i);
             b(i)=-D(i)+D_w(i);
             c(i)=-D(i);
        end
    end        

    b(end)=[];
    b(1)=b(1)-D(1);
    c(1)=[];
    y=T_initial-((rho_i*L./((phi_np1_km1_j*(rho_br*c_br)+...
        (1-phi_np1_km1_j)*(rho_i*c_i)))).*...
        (phi_np1_km1_j-phi_initial))+(1/((phi_np1_km1_j*(rho_br*c_br)+...
        (1-phi_np1_km1_j)*(rho_i*c_i))))'.*T_source;
    if w_np1_k_j<0
        y(end)=y(end)+D(end)*Tbottom-D_w(end)*Tbottom;
    else
        y(end)=y(end)+D(end)*Tbottom;
    end
    if BC==2;
        y(1)=y(1)-2*D(1)*dz*T_grad;
    else
        y(1)=Ttop;
    end
    % Solve for T
    T_np1_k_j=Thomas_Trid(a,b,c,y');
    
    T_np1_km1_j=T_np1_k_j;
    
    %% Create K_s matrix and solve for S
    one_over_phi1(1:matrix_dimension)=1;
    one_over_phi=rdivide(one_over_phi1',phi_np1_km1_j);
    for i=1:matrix_dimension;
        if one_over_phi(i)==Inf;
            one_over_phi(i)=0;
        end;
    end;
    k_s_mod1(1:matrix_dimension)=k_s;
    k_s_mod=times(k_s_mod1',phi_np1_km1_j.^m);
    D_s_mod=((k_s_mod*dt)/dz^2).*one_over_phi;

    D_w_s=((w_np1_k_j*dt)/(dz)).*one_over_phi;

    a2(1)=1+D_s_mod(1);
    b2(1)=-D_s_mod(1);
    for i=2:matrix_dimension;
        if w_np1_k_j(i)==0;
            a2(i)=1+2*D_s_mod(i);
            b2(i)=-D_s_mod(i);
            c2(i)=-D_s_mod(i);
        elseif w_np1_k_j(i)>0;
            a2(i)=1+2*D_s_mod(i)+D_w_s(i);
            b2(i)=-D_s_mod(i);
            c2(i)=-D_s_mod(i)-D_w_s(i);
        else
            a2(i)=1+2*D_s_mod(i)-D_w_s(i);
            b2(i)=-D_s_mod(i)+D_w_s(i);
            c2(i)=-D_s_mod(i);
        end
    end

    b2(end)=[];
    c2(1)=[];
    y2=S_initial-(rho_i/rho_br)*one_over_phi.*S_initial...
        .*(phi_np1_km1_j-phi_initial)+one_over_phi.*S_source;
    if w_np1_k_j(end)<0
        y2(end)=y2(end)+D_s_mod(end)*S_bottom+D_w_s(end)*S_bottom;
        y2(1)=y2(1);
    else
        y2(end)=y2(end)+D_s_mod(end)*S_bottom;
        y2(1)=y2(1);
    end

    %Solve for S
    S_np1_k_j=Thomas_Trid(a2,b2,c2,y2');
    
        S_np1_km1_j=S_np1_k_j;
    
    %% Appending value to matrix to check for convergence
    S_evolve=[S_evolve S_np1_k_j];
    T_evolve=[T_evolve T_np1_k_j];
    phi_evolve=[phi_evolve phi_np1_km1_j];
    
    %% Reassign T vector to use in next k iteration
    T_np1_km1_j=T_np1_k_j;
    phi_np1_k_j=phi_np1_km1_j;
    % tolerance tests
    if counter==1;
        SErr=SErr;
        TErr=TErr;
        PhiErr=PhiErr;
    else
        TErr=max(abs(T_evolve(:,counter)-T_evolve(:,counter-1)));
        PhiErr=max(abs(phi_evolve(:,counter)-phi_evolve(:,counter-1)));
        SErr=max(abs(S_evolve(:,counter)-S_evolve(:,counter-1)));
    end;
end



