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
%Tol=.001
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
    
    %delta_T=1.853*(S_np1_km1_j/28);        %% Old Freezing point depression due to salt
    delta_T=Tm-(-(1.333489497*(10^-5)*S_np1_km1_j.^2)-0.01612951864*S_np1_km1_j+...
        273.055175687);                       %% Liquidus point using FREZCHEM Europa ocean
%     delta_T=Tm-(-(9.1969758*(10^-5)*S_np1_km1_j.^2)-0.03942059*S_np1_km1_j+...
%        272.63617665);                        %% Liquidus point using FREZCHEM seawater
    Hs=c_i*(Tm-delta_T);                     %% Enthalpy Calculation
    Hsmelt=c_i*Tm;                           %% Enthalpy to melt ice
    %% Calculating solid fraction at k-1 iteration for use in
    %% temperature, salinity, and brine velocity at k iteration
    
    %% Value of Phi(n+1,k-1,j)
    %phi_np1_km1_j=((c_br*T_np1_km1_j+phi_np1_km2_j*L)-Hs)./...
    %    (c_br*(Tm-delta_T)-c_i*(Tm-delta_T));
    phi_np1_km1_j=phi_np1_km2_j;
    for i=1:matrix_dimension;
    %% Test to see if melting ice or freezing brine
    %if T_np1_km1_j(i)<Tm-delta_T(i);
    En=c_i*T_np1_km1_j(i)+L*phi_np1_km2_j(i);
        if En<Hs(i)
            %phi_np1_km1_j(i)<0;
            phi_np1_km1_j(i)=0;
        elseif En>Hs(i)+L;
            phi_np1_km1_j(i)=1;%phi_np1_km2_j(i);
        else
            phi_np1_km1_j(i)=((c_i*T_np1_km1_j(i)+phi_np1_km2_j(i)*L)-Hs(i))/...
            L;
            %(c_br*(Tm-delta_T(i))-c_i*(Tm-delta_T(i)));
            %phi_np1_km1_j(i)=phi_np1_km1_j(i);
        end;
%     else
%     En=c_i*T_np1_km1_j(i)+L*phi_np1_km2_j(i);
%         if En>Hsmelt+L;
%             phi_np1_km1_j(i)=1;
%         else
%             phi_np1_km1_j(i)=((c_i*T_np1_km1_j(i)+phi_np1_km2_j(i)*L)-Hsmelt)/...
%             L;
%         end;
%     end;
    end;
    
    %% limiting the maximum freeze out to avoid over saliniation
    
    %mark=0*[0:dz:H];
    %for i=1:matrix_dimension
    %    if phi_initial(i)<=0.03 && phi_np1_km1_j(i)<phi_initial(i);
    %        phi_np1_km1_j(i)=phi_initial(i);
    %        mark(i)=1;
    %    %elseif phi_np1_km1_j(i)==0;
    %    %    phi_np1_km1_j(i)=phi_initial(i);
    %    else
    %    end
    %end
    
    %cycle=matrix_dimension;
    %for i=matrix_dimension:-1:1;
    %    if phi_initial(i)>=0.02;
    %        cycle=cycle-1;
    %    else
    %    continue
    %    end
    %end
    %if cycle==1;
    %    continue
    %else
    %        for j=1:cycle;
    %            if phi_np1_km1_j(j)<phi_initial(j);
    %                phi_np1_km1_j(j)=phi_initial(j);
    %            else
    %            end
    %        end
    %end
    
    
    %if freeze==1;
    %for i=1:matrix_dimension;
    %    if phi_np1_km1_j(i)<min(phi_initial);
    %        phi_np1_km1_j(i)=phi_initial(i);
    %    else
    %        phi_np1_km1_j(i)=phi_np1_km1_j(i);
    %    end
    %end
    %else
    %end
    
        
    
    %% Freeze out of elements with liquid frac. < 5%
    %if n>0;
    %for i=1:matrix_dimension;
    %    if phi_np1_km1_j(i)<.05;
    %        phi_np1_km1_j(i)=0;
    %    else
    %        phi_np1_km1_j(i)=phi_np1_km1_j(i);
    %    end
    %end
    %else
        phi_np1_km1_j=phi_np1_km1_j;
    %end
    %% Reassigning Phi(n+1,k-2) as Phi(n+1,k-1) for next
    %% iteration
    phi_np1_km2_j=phi_np1_km1_j;
    
    %% Solving for Darcy Velocity (eq. 9 Hunke et. al.) no gravity drainage
    %aw(1:matrix_dimension)=1;
    %bw(1:matrix_dimension-1)=-1/2;
    %cw(1:matrix_dimension-1)=-1/2;
    
%     if n<2;
         w_np1_k_j=w_initial;
%     else
%     for i=2:matrix_dimension;
%          if phi_np1_km1_j(i)<=0.05;
%             w_np1_k_j_loop(1:i)=0;
%          else
%             w_np1_k_j_loop(i)=w_np1_k_j_loop(i-1)+(dz/dt)*(rho_i/rho_br-1)*...
%                 (phi_np1_km1_j(i)-phi_initial(i));
%          end;
%     end;
    
    %count=matrix_dimension;
    %for i=matrix_dimension:-1:2;
    %    if abs(w_np1_k_j_loop(i))>0;
    %        count=count-1;
    %    else
    %        break
    %    end
    %end
    %w_np1_k_j_loop(1:count)=0;
    
    %(((3*10^-8)*phi_np1_km1_j.^2)./(1.787*10^-3))...
    %.*((50/68)*(S_initial-34)*9.8)

    w_np1_k_j=w_np1_k_j_loop;
    w_np1_k_j_loop=w_initial;
    
    %end

    %w_np1_k_j=Thomas_Trid(aw,bw,cw,yw')
    
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
    %a=1+2*D;
    %b=(D_w-D);
    %c=-D_w-D;
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
    %for i=1:matrix_dimension;
    %    if mark(i)==1;
    %        D_s_mod(i)=0;
    %    else
    %    end
    %end
    D_w_s=((w_np1_k_j*dt)/(dz)).*one_over_phi;
    %a2=1+2*D_s_mod.*one_over_phi;
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
    %b2=-one_over_phi.*(D_w_s-D_s_mod);
    %c2=one_over_phi.*(-D_w_s-D_s_mod);
    b2(end)=[];
    c2(1)=[];
    y2=S_initial-(rho_i/rho_br)*one_over_phi.*S_initial...
        .*(phi_np1_km1_j-phi_initial)+one_over_phi.*S_source;
    if w_np1_k_j(end)<0
        y2(end)=y2(end)+D_s_mod(end)*S_bottom+D_w_s(end)*S_bottom;
        y2(1)=y2(1);%+D_s_mod(1)*S_np1_km1_j(1);
    else
        y2(end)=y2(end)+D_s_mod(end)*S_bottom;
        y2(1)=y2(1);%+D_s_mod(1)*S_np1_km1_j(1);
    end
    
    % Solve for S
    %S_np1_k_j=Thomas_Trid(a2,b2,c2,y2')
    
    % Modifying the S tensor so it isn't singular (i.e. neglecting zero
    % filled rows)
%    point=1;
%    for i=1:matrix_dimension;
%        if a2(i)==0;
%            point=point+1;
%        else point=point;
%        end
%    end
%    mod_S=zeros(matrix_dimension);
%    for i=1:matrix_dimension;
%        mod_S(i,i)=a2(i);
%    end
%    for i=1:matrix_dimension-1;
%        mod_S(i,i+1)=b2(i);
%        mod_S(i+1,i)=c2(i);
%    end
%    mod_S2=mod_S(point:matrix_dimension,point:matrix_dimension);
%    Smat=zeros(matrix_dimension);
%    Smat(point:matrix_dimension,point:matrix_dimension)=inv(mod_S2);
    % above line needs cleaned up to not use inv() command
    
%    S_np1_k_j=Smat*y2;
    S_np1_k_j=Thomas_Trid(a2,b2,c2,y2');
    
    %if n==1
    %    S_np1_k_j=S_initial;
    %else
    %    S_np1_k_j=S_np1_k_j;
    %end
    
    % locating the element that salt from frozen out cells will be sent to
%    solid_point=0;
%    salts=0;
%    for i=1:matrix_dimension;
%        if S_np1_k_j(i)==0;
%            solid_point=solid_point+1;
%            salts=salts+phi_initial(i)*S_initial(i);
%        end;
%    end;
%    expulsion_element=matrix_dimension;
%    for i=matrix_dimension:-1:1;
%        if phi_np1_km1_j(i)>0;
%            expulsion_element=expulsion_element-1;
%        else
%            break
%        end
%    end
    %S_np1_k_j(1:expulsion_element)=0
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
        %Err=max([Err_T,Err_phi,Err_S]);
    end;
% brine expulsion for next iteration
%S_np1_k_j(expulsion_element+1)=S_np1_k_j(expulsion_element+1)+phi_np1_km1_j(expulsion_element+1)*salts;
end
%S_np1_k_j(expulsion_element+1)=S_np1_k_j(expulsion_element+1)+phi_np1_km1_j(expulsion_element+1)*salts;



