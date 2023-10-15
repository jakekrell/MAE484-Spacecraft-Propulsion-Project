%% Development Information
% MAE 484 Spacecraft Propulsion 
% CreateNozzle2D_G10.m
% 
% Create struct with information regarding 2D drawing of nozzle.
%
% input[]: description. type
%
% output[Nozzle]: information regarding 2D drawing of nozzle. struct
% 
% Assumptions:
% (1) 1.5*rt   == chamber to conv. section, and conv. to throat inlet
% (2) 0.382*rt == throat outlet to exit section
% (3) origin at center of throat
% 
% 
% Primary Developer Contact Information:
% Jacob P. Krell [Project Group 10]
% Aerospace Engineering Undergraduate Student
% Statler College of Engineering & Mineral Resources
% Dept. Mechanical and Aerospace Engineering
% West Virginia University (WVU)
% jpk0024@mix.wvu.edu
%
%
%
% Development History
% Date              Developer        Comments
% ---------------   -------------    --------------------------------
% Oct. 14, 2023     J. Krell         For Spacecraft Propulsion Project
%
%%

function Nozzle = CreateNozzle2D_G10(rc,Lc,rt1,rt,rt2,a,L,res,plt,xc_offset)

    if (rc-rt) > 2*rt1 % chamber to throat cannot be connected smoothly
        error('Cannot connect geometry. Increase rt and/or rt1, and/or decrease rc.')
        elseif rc < rt
        error('Chamber radius must be greater than throat radius.')
    end

    a = a*pi/180; % convert alpha to radians
    x = linspace(0,1,res); % x-grid spacing per curve

    % find where throat ends
    
%     % equations:
%     y2 = rt + rt2 - sqrt(rt2^2 - x2^2);
%     yL = tan(a) * (x - x2_end) + y2_end;
% 
%     % need derivatives equal at (x2_end,y2_end):
%     dy2 = -1/2 * (rt2^2 - x2^2)^(-1/2) * (-2*x2);
%     dyL = tan(a);
%
%     % solving for x2_end:
%     dy2 = dyL --> x2_end/sqrt(rt2^2-x2_end^2) = tan(a);
    x2_end = rt2 * tan(a) / sqrt(1 + tan(a)^2);
    
    x2 = x*x2_end;
    y2 = rt + rt2 - sqrt(rt2^2 - x2.^2);

    xL = [x2_end, L];
    yL = [y2(end), tan(a)*(L-x2_end)+y2(end)];


    % find where throat starts:

%     % equations:
%     y1 = rt + rt1 - sqrt(rt1^2 - x1^2);
%     y0 = sqrt(rc^2 - (x0-xc_end)^2);
% 
%     % need derivatives equal at (x1_end,y1_end):
%     dy1 = -1/2 * (rt1^2 - x1^2)^(-1/2) * (-2*x1);
%     dy0 = 1/2 * (rc^2 - (x0-xc_end)^2))^(-1/2) * (-2*(x0-xc_end));
%
%     % solving for xc_end given dy1 = dy0 at x1_end = x0_end:
%     dy1 = dy0 -->
%     -x1/sqrt(rt1^2-x1^2) = (x0-xc_end)/sqrt(rc^2-(x0-xc_end)^2)); -->
%     -x1_end/sqrt(rt1^2-x1_end^2) . . .
%     . . . = (x0_end-xc_end)/sqrt(rc^2-(x0_end-xc_end)^2); ......... (eq1)

%     % solving for xc_end given y1 = y0 at x1_end = x0_end:
%     y1 = y0 -->
%     rt+rt1-sqrt(rt1^2-x1^2) = sqrt(rc^2-(x0-xc_end)^2); -->
%     rt+rt1-sqrt(rt1^2-x1_end^2) = sqrt(rc^2-(x0_end-xc_end)^2); ... (eq2)
    syms x1_end xc_end
    eq1 = -x1_end/sqrt(rt1^2-x1_end^2) == (x1_end-xc_end)/sqrt(rc^2-(x1_end-xc_end)^2);
    eq2 = rt+rt1-sqrt(rt1^2-x1_end^2) == sqrt(rc^2-(x1_end-xc_end)^2);
    
    ig_y = (rc+rt)/2;
    ig_x1 = -sqrt(rt1^2-ig_y^2);
    ig_xc = ig_x1 - sqrt(rc^2-ig_y^2);
    sol = vpasolve([eq1,eq2],[x1_end,xc_end],[ig_x1;ig_xc]);
    
    x1_end = double(sol.x1_end);
    xc_end = double(sol.xc_end);

    if x1_end > 0
        x1_end = -x1_end;
    end
    if xc_end > 0
        xc_end = -xc_end;
    end
    size_x1_end = size(x1_end); size_xc_end = size(xc_end);
    if size_x1_end(1) == 0 || size_xc_end(1) == 0
        error('A numerical solution could not be found for the chamber to contraction location.')
    end
    
    x1 = flip(x)*x1_end;
    y1 = rt + rt1 - sqrt(rt1^2 - x1.^2);

    x0 = x*(x1_end-xc_end) + xc_end;
    y0 = sqrt(rc^2 - (x0-xc_end).^2);

    xc = [xc_end-Lc, xc_end];
    yc = [rc, rc];


    % combine sections

    x_all = [xc,x0(2:end),x1(2:end),x2(2:end-1),xL];
    y_all = [yc,y0(2:end),y1(2:end),y2(2:end-1),yL];

    if xc_offset
        offset = xc(1);
        x_all = x_all - offset;
        xc = xc - offset;
        x0 = x0 - offset;
        x1 = x1 - offset;
        x2 = x2 - offset;
        xL = xL - offset;
    end

    Nozzle.x = x_all';
    Nozzle.y = [-y_all;y_all]';
    Nozzle.Chamber.x = xc';
    Nozzle.Chamber.y = [-yc;yc]';
    Nozzle.Contraction.x = x0';
    Nozzle.Contraction.y = [-y0;y0]';
    Nozzle.Throat.Upstream.x = x1';
    Nozzle.Throat.Upstream.y = [-y1;y1]';
    Nozzle.Throat.Downstream.x = x2';
    Nozzle.Throat.Downstream.y = [-y2;y2]';
    Nozzle.Expansion.x = xL';
    Nozzle.Expansion.y = [-yL;yL]';

    Nozzle.Characteristics.Lc = Lc;
    Nozzle.Characteristics.L = L;
    Nozzle.Characteristics.alpha = a;
    Nozzle.Characteristics.rc = rc;
    Nozzle.Characteristics.rt1 = rt1;
    Nozzle.Characteristics.rt = rt;
    Nozzle.Characteristics.rt2 = rt2;
   

    % plotting:

    if plt{1}
        PlotNozzle2D_G10(Nozzle,plt)
    end

end


function PlotNozzle2D_G10(Nozzle,plt)
    
    rt1 = Nozzle.Characteristics.rt1;
    rt = Nozzle.Characteristics.rt;
    rt2 = Nozzle.Characteristics.rt2;
    
    x_all = Nozzle.x;
    y_all = Nozzle.y(:,2);
    xc = Nozzle.Chamber.x;
    yc = Nozzle.Chamber.y(:,2);
    x0 = Nozzle.Contraction.x;
    y0 = Nozzle.Contraction.y(:,2);
    x1 = Nozzle.Throat.Upstream.x;
    y1 = Nozzle.Throat.Upstream.y(:,2);
    x2 = Nozzle.Throat.Downstream.x;
    y2 = Nozzle.Throat.Downstream.y(:,2);
    xL = Nozzle.Expansion.x;
    yL = Nozzle.Expansion.y(:,2);

    plt_global = plt{1};
    plt_sect = plt{2};
    plt_mid = plt{3};
    plt_radii = plt{4};
    plt_nTyp = plt{5};
    plt_nPnt = plt{6};
    plt_sTyp = plt{7};
    plt_sPnt = plt{8};
    plt_mTyp = plt{9};
    plt_mPnt = plt{10};
    plt_rTyp = plt{11};
    plt_rPnt = plt{12};
    plt_buffer = plt{13}; % percentage of 1 side's negative space vs image

    if plt_global
        f1 = figure(1);
        axis equal
        hold on
        grid on

        plot(x_all,y_all,plt_nTyp,LineWidth=plt_nPnt)
        plot(x_all,-y_all,plt_nTyp,LineWidth=plt_nPnt)

        if plt_sect
            plot([xc(1),xc(1)],[-yc(1),yc(1)],plt_sTyp,LineWidth=plt_sPnt)
            plot([x0(1),x0(1)],[-y0(1),y0(1)],plt_sTyp,LineWidth=plt_sPnt)
            plot([x1(1),x1(1)],[-y1(1),y1(1)],plt_sTyp,LineWidth=plt_sPnt)
            plot([x2(1),x2(1)],[-y2(1),y2(1)],plt_sTyp,LineWidth=plt_sPnt)
            plot([xL(1),xL(1)],[-yL(1),yL(1)],plt_sTyp,LineWidth=plt_sPnt)
            plot([xL(end),xL(end)],[-yL(end),yL(end)],plt_sTyp,LineWidth=plt_sPnt)
        end

        if plt_mid
            plot([xc(1),xL(end)],[0,0],plt_mTyp,LineWidth=plt_mPnt)
        end

        if plt_radii
            for pm = [1,-1]
                plot([x0(1),x0(1)],pm*[0,y0(1)],plt_rTyp,LineWidth=plt_rPnt)
                plot([x0(1),x0(end)],pm*[0,y0(end)],plt_rTyp,LineWidth=plt_rPnt)
                plot([x1(end),x1(1)],pm*[rt+rt1,y1(1)],plt_rTyp,LineWidth=plt_rPnt)
                plot([x1(end),x1(end)],pm*[rt+rt1,y1(end)],plt_rTyp,LineWidth=plt_rPnt)
                plot([x2(1),x2(1)],pm*[rt+rt2,y2(1)],plt_rTyp,LineWidth=plt_rPnt)
                plot([x2(1),x2(end)],pm*[rt+rt2,y2(end)],plt_rTyp,LineWidth=plt_rPnt)
            end
        end

        maxY = max(y_all);
        if plt_buffer > 0
            bufferX = plt_buffer*(x_all(end)-x_all(1))/(1-2*plt_buffer);
            bufferY = 2*plt_buffer*maxY/(1-2*plt_buffer);
            xlim([x_all(1)-bufferX, x_all(end)+bufferX])
            ylim([-maxY-bufferY, maxY+bufferY])
        else
            xlim([x_all(1), x_all(end)])
            ylim([-maxY, maxY])
        end

        title('Conical Rocket Nozzle Design','Interpreter','latex')
        xlabel('Location Along Nozzle Length (m)','Interpreter','latex')
        ylabel('Location Along Nozzle Height (m)','Interpreter','latex')
        hold off
        set(gca,'TickLabelInterpreter','latex')

    end
end

