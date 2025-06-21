% Copyright 2021 Delft University of Technology
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%
%      http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

function Visualisation(t,p_t_x,p_t_y,p_t_z, pos_O,eta_OB, Params, sampletime, view_az, view_el, subplot_num)   
%     persistent handleTether
%     persistent pathpoints
%     persistent Vert
%     persistent handlePlane
    persistent fig_prop
    global bInitConditionflag
    windDirection_rad = Params.env.wind.direction;       
%     coder.extrinsic('plot3','patch','animatedline','addpoints');
%     t = others.time(k);
    subplot(2,2,subplot_num)
    %% Figure settings
    if t==0
        
%           addToolbarExplorationButtons(gcf)
        xlabel('x_W (m)');%, 'interpreter', 'latex')
        ylabel('y_W (m)');%, 'interpreter', 'latex')
        zlabel('z_W (m)');%, 'interpreter', 'latex')

        grid on
        box on;
        axis equal; hold on;

%         view(30,21); % set the view angle for figure
         view(view_az,view_el)
         ax = [-50 800];
         ay = [-450 450];
         az = [0 700];
         axis([ax,ay,az]);  hold on

         X1 = [ax(1);ax(2);ax(2);ax(1)];
         Y1 = [ay(1);ay(1);ay(2);ay(2)];
         Z1 = [0;0;0;0];

         c = [0 1 0; % red
             0 1 0; % green
             0 1 0;
             0 1 0];
         fill3(X1,Y1,Z1, c(:,2),'EdgeColor','none','FaceColor',c(1,:),'FaceAlpha',0.1);
         fig_prop(subplot_num).handleTether = [];
         fig_prop(subplot_num).pathpoints = [];
         fig_prop(subplot_num).handlePlane = [];
    end

    %% Draw tether particles
    if isempty(fig_prop(subplot_num).handleTether)

        M_WO = [cos(windDirection_rad), sin(windDirection_rad), 0;
            sin(windDirection_rad), -cos(windDirection_rad), 0;
            0, 0, -1];
        pos_W = M_WO*[pos_O(1);pos_O(2);pos_O(3)];

%         p_t_x = X.T(1:3:45,k);
%         p_t_y = X.T(2:3:45,k); 
%         p_t_z = X.T(3:3:45,k);
        handleTether = plot3( [0;flipud(p_t_x);pos_W(1)], [0;flipud(p_t_y);pos_W(2)],[0;flipud(p_t_z);pos_W(3)],...
            '-', 'color', [0.1 0.1 0.1], 'Markersize',.02, 'Linewidth', .02,...
            'MarkerFaceColor', [0.1 0.1 0.1]); hold on;
        fig_prop(subplot_num).handleTether = handleTether;
    end

    %% Draw aircraft
    if isempty(fig_prop(subplot_num).handlePlane)
        pn = pos_O(1);
        pe = pos_O(2);
        pd = pos_O(3);
        phi = eta_OB(1);
        theta = eta_OB(2);
        psi = eta_OB(3);

        M_WO = [cos(windDirection_rad), sin(windDirection_rad), 0;
            sin(windDirection_rad), -cos(windDirection_rad), 0;
            0, 0, -1];

        scale = 1.5;
        %         [F, Vert, C] = rndread('plane_thesis_big.stl');
        Image_kite = coder.load('plane_image.mat','F','Verti','C');
        F = Image_kite.F;
        Verti = Image_kite.Verti;
        C = Image_kite.C;
        fig_prop(subplot_num).Vert = scale*Verti';
        V = [ cos(psi)*cos(theta) , cos(psi)*sin(phi)*sin(theta) - cos(phi)*sin(psi) , sin(phi)*sin(psi) + cos(phi)*cos(psi)*sin(theta);
            cos(theta)*sin(psi)   , cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(theta) , cos(phi)*sin(psi)*sin(theta) - cos(psi)*sin(phi);
            -sin(theta)        ,             cos(theta)*sin(phi)                  ,             cos(phi)*cos(theta)                 ]...
            * fig_prop(subplot_num).Vert;
        V = V + [pn*ones(1,size(V,2));pe*ones(1,size(V,2));pd*ones(1,size(V,2))];
        V = M_WO*V;
%         
        if bInitConditionflag == 1 
            col =  [ 1 0, 0  ];
        else 
            col =  [ 57/255, 106/255, 177/255  ];
        end
        handlePlane = patch('faces', F, 'vertices' ,V');
        set(handlePlane, 'facec', 'r');              % Set the face color (force it)
        set(handlePlane, 'FaceVertexCData', C);       % Set the color (from file)
        set(handlePlane, 'EdgeColor','none');
        set(handlePlane, 'FaceLighting', 'none');
        set(handlePlane, 'FaceColor', col );
        fig_prop(subplot_num).handlePlane = handlePlane;
    end

    %% Draw black thin trail behind kite
    if isempty(fig_prop(subplot_num).pathpoints)
        if bInitConditionflag == 1 
            pathpoints = animatedline('Linestyle','-','color', [0 0 0], 'Linewidth', 1);
            fig_prop(subplot_num).pathpoints = pathpoints;
            
        else
            pathpoints = animatedline('Linestyle','--','color', [0 0 0], 'Linewidth', 1);
            fig_prop(subplot_num).pathpoints = pathpoints;
        end
    end



    if (mod(t,sampletime)<=10e-2 && t~=0)
        %% Change each element on figure to current positions
        pn = pos_O(1);
        pe = pos_O(2);
        pd = pos_O(3);
        phi = eta_OB(1);
        theta = eta_OB(2);
        psi = eta_OB(3);

        M_WO = [cos(windDirection_rad), sin(windDirection_rad), 0;
            sin(windDirection_rad), -cos(windDirection_rad), 0;
            0, 0, -1];

        pos_W = M_WO*[pos_O(1);pos_O(2);pos_O(3)];

        % Tether
%         p_t_x = X.T(1:3:45,k);
%         p_t_y = X.T(2:3:45,k); 
%         p_t_z = X.T(3:3:45,k);
        set(fig_prop(subplot_num).handleTether,'XData',[0;flipud(p_t_x);pos_W(1)], 'YData',[0;flipud(p_t_y);pos_W(2)],'ZData',[0;flipud(p_t_z);pos_W(3)]);

        % Kite
        V = [ cos(psi)*cos(theta), cos(psi)*sin(phi)*sin(theta) - cos(phi)*sin(psi), sin(phi)*sin(psi) + cos(phi)*cos(psi)*sin(theta);
            cos(theta)*sin(psi), cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(theta), cos(phi)*sin(psi)*sin(theta) - cos(psi)*sin(phi);
            -sin(theta),                              cos(theta)*sin(phi),                              cos(phi)*cos(theta)] * fig_prop(subplot_num).Vert;
        V = V + [pn*ones(1,size(V,2));pe*ones(1,size(V,2));pd*ones(1,size(V,2))];
        V = M_WO*V;
        set(fig_prop(subplot_num).handlePlane, 'Vertices', V');
        % Trail
        pos_W = M_WO*[pn;pe;pd];
        addpoints(fig_prop(subplot_num).pathpoints,pos_W(1), pos_W(2),pos_W(3))
        drawnow;
        str = sprintf('simulation time %.2f sec', t);
        title (str)
    end
 end