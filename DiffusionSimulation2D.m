function [DiffOut1,DiffOut2,DiffOut3,DiffOut4,DiffOut5,DiffOut6,DiffOut7,DiffOut8,DiffOut9] = DiffusionSimulation2D(H2O2_Conc,OPD_Conc,Diff_Desc)
    if nargin ~= 3
        disp('wrong number of inputs');
        disp('such as (0.5,0.5,[1 1]) as input');
        return
    end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%setup field
% well size = 2cm x 2cm.
% matrix of 364 x 364 was used.
% solution volume = 1 ml.
% each grid = 55um x 5um corresponding to average size of coacevate.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Define H2O2 field
    if Diff_Desc(1) == 1 % diffuse H2O2 from one edge of field
        field1 = zeros(364,364); 
        for i = 1:length(field1)
            for j = 1:length(field1)
                if ((i-364)^2 + (j-182)^2)*3025 <= 1.6e7/pi 
                    field1(i,j) = 50*H2O2_Conc*0.5/0.5053;%mM H2O2
                end
            end
        end
    elseif Diff_Desc(1) == 0 % homogenous field 
        field1 = OPD_Conc*ones(364,364);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%Define o-PD field
    if Diff_Desc(2) == 1
        field2 = zeros(364,364); % diffuse o-PD from one edge of field
        for i = 1:length(field2)
            for j = 1:length(field2)
                if ((i-364)^2 + (j-182)^2)*3025 <= 1.6e7/pi 
                    field2(i,j) =50*OPD_Conc*0.5/0.5053;%mM aminophenol
                end
            end
        end
    elseif Diff_Desc(2) == 0 % homogenous field 
        field2 = OPD_Conc*ones(364,364);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define Coacevate field
    %size of coacevate array. Here is 15 x 15.
    Coa_Num = 15;
    %HPR enzyme concentration in coacevate
    HPR = 3.66e-7.*ones(Coa_Num,Coa_Num); 
    % Construct intermediates' field
    HPR_H2O2 = zeros(Coa_Num,Coa_Num);
    CompoundI = zeros(Coa_Num,Coa_Num);
    CompoundI_OPD = zeros(Coa_Num,Coa_Num);
    CompoundII = zeros(Coa_Num,Coa_Num);
    CompoundII_OPD = zeros(Coa_Num,Coa_Num);
    DAP = zeros(Coa_Num,Coa_Num);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
    %setup diffusion parameters
    nodewidth = 5.5e-5; %m(meter)
    D1 = 1.71e-9;%m^2/s for H2O2;
    D2= 6.57e-10; %m^2/s for o-PD
    dt = 0.1;%s. step of simulation
    % in total of 3h
    Tend = 180*60/dt;
    % Compute Fo factor
    instability = D1*dt/(nodewidth)^2;
    if instability >=0.25
        disp('Fo = ',num2str(instability));
        disp('consider to reduce simulation step for stable simulation');
        % key diffusion factor, above 0.25 at 2D will lead to instability of
        % system
        return
    end
    % recording step for concentration change in regard to diffusion direction.
    j = 1:8/dt:(Tend+10/dt); 
    % recording step for whole chemical field.
    j2 = 1:60/dt:(Tend+10/dt);
    % initialization of data collection
    k = 1;
    k2 = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Chemical diffusion + enzymatic conversion
    for i = 1:Tend
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %record substrate concentrations regarding diffusion direction.
         if isequal(j(k),i)
           time(k) = (j(k)-1)*dt;
           Conc_record_1_row(:,k) = field1(182,[182-Coa_Num+1:2:182+Coa_Num+1]);
           Conc_record_2_row(:,k) = field2(182,[182-Coa_Num+1:2:182+Coa_Num+1]);
           Conc_record_1_col(:,k) = field1([182-Coa_Num+1:2:182+Coa_Num+1],182);
           Conc_record_2_col(:,k) = field2([182-Coa_Num+1:2:182+Coa_Num+1],182);         
           DAP_record(:,:,k) = DAP;
           k = k+1;
         end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %record whole chemical field.
        if isequal(j2(k2),i)
           %time2(k2) = (j2(k2)-1)*dt;
           H2O2_field(:,:,k2) = field1;
           OPD_field(:,:,k2) = field2;
           k2 = k2+1 % to notice the progress of simulation
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Diffusion equations
        % compute chemical 1 (H2O2) diffusion
        f1_yminush = [field1(2:end,:);field1(end,:)];
        f1_xminush = [field1(:,2:end) field1(:,end)];

        f1_yplush = [field1(1,:);field1(1:end-1,:)];
        f1_xplush = [field1(:,1) field1(:,1:end-1)];

        d_f1 = (f1_yplush + f1_xplush + f1_yminush + ...
                f1_xminush - 4*field1)./(nodewidth^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute chemical 2 (o-PD) diffusion
        f2_yminush = [field2(2:end,:);field2(end,:)];
        f2_xminush = [field2(:,2:end) field2(:,end)];

        f2_yplush = [field2(1,:);field2(1:end-1,:)];
        f2_xplush = [field2(:,1) field2(:,1:end-1)];

        d_f2 = (f2_yplush + f2_xplush + f2_yminush + ...
                f2_xminush - 4*field2)./(nodewidth^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do one diffusion step
        d1 = D1*d_f1*dt;
        d2 = D2*d_f2*dt;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%% enzymatic reactions in quasi-equilibrium fashion
% do reaction on every defined coacevate in the diffusion step time 
         for m = 1:Coa_Num
               for n = 1:Coa_Num
                   [HPR(m,n),field1(180-Coa_Num+1+2*m,180-Coa_Num+1+2*n),...
                       HPR_H2O2(m,n),CompoundI(m,n),field2(180-Coa_Num+1+2*m,...
                       180-Coa_Num+1+2*n),CompoundI_OPD(m,n),CompoundII(m,n)...
                       ,CompoundII_OPD(m,n),DAP(m,n)] = ...
                       HRP_pingpong_Complicated(HPR(m,n),...
                       field1(180-Coa_Num+1+2*m,180-Coa_Num+1+2*n),...
                       HPR_H2O2(m,n),CompoundI(m,n),...
                       field2(180-Coa_Num+1+2*m,180-Coa_Num+1+2*n),...
                       CompoundI_OPD(m,n),CompoundII(m,n),...
                       CompoundII_OPD(m,n),DAP(m,n),dt);
               end
         end
% update the change to chemical field
        field1 = field1+ d1;
        field2 = field2+ d2;
% feed to next round of simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%x`%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Construct outputs
    DiffOut1 = Conc_record_1_col; % H2O2 in column
    DiffOut2 = Conc_record_1_row; % H2O2 in row
    DiffOut3 = Conc_record_2_col; % o-PD in column
    DiffOut4 = Conc_record_2_row; % o-PD in row
    DiffOut5 = H2O2_field;        % whole H2O2 diffusion field
    DiffOut6 = OPD_field;         % whole o-PD diffusion field
    DiffOut7 = j2;                % time step of data recording 
    DiffOut8 = j;                 % time step of data recording
    DiffOut9 = DAP_record;        % final concentration of radicals in all coacevates
end