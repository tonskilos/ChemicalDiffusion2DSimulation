
function [OUT1,OUT2,OUT3,OUT4,OUT5,OUT6,OUT7,OUT8,OUT9] = HRP_pingpong_Complicated(V1,V2,V3,V4,V5,V6,V7,V8,V9,t)
    initial = [V1,V2,V3,V4,V5,V6,V7,V8,V9];
    k1 = 2.7e3;%mM-1s-1
    k2 = 2170; % s-1
    k3 = 7.1e5; % mM-1s-1
    k4 = 2.17e5; % s-1
    k5 = 7.1e3;%mM-1s-1
    k6 = 437; % s-1
    sol = ode23s(@func,[0 t],initial);
    OUT1 = sol.y(1,end); OUT2 = sol.y(2,end);
    OUT3 = sol.y(3,end); OUT4 = sol.y(4,end);
    OUT5 = sol.y(5,end); OUT6 = sol.y(6,end);
    OUT7 = sol.y(7,end); OUT8 = sol.y(8,end);
    OUT9 = sol.y(9,end);
    %plot(sol.x,sol.y);
    %legend('HPR','H2O2','HPR_oxi','OPD','DAP')
    function dydt = func(t,y)
        dydt = zeros(9,1);
        dydt(1) = -k1*y(1)*y(2)+k6*y(8);%HPR
        dydt(2) = -k1*y(1)*y(2)-dydt(9);%H2O2
        dydt(3) = k1*y(1)*y(2) - k2*y(3); %HPR-H202
        dydt(4) = k2*y(3) - k3*y(4)*y(5); % Compound I
        dydt(5) = -k3*y(4)*y(5) - k5*y(5)*y(7); %OPD
        dydt(6) = k3*y(4)*y(5) - k4*y(6); % Compound I-OPD
        dydt(7) = k4*y(6) - k5*y(5)*y(7); % Compound II
        dydt(8) = k5*y(5)*y(7) - k6*y(8); % Compound II-OPD
        dydt(9) = k6*y(8)+ k4*y(6); %Radicals
    end
end