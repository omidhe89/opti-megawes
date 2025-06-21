function animation_awe = awe_animation(t, optimal_results, h, params_num)
    
 % Copyright (c) 2025 Omid Heydarnia
%
% This software is licensed under the Apache License, Version 2.0 (the "License");
% See the accompanying LICENSE file for the full text of the license.
%
% Author: Omid Heydarnia
% Email: omid.heydarnia@ugent.be
% Date: June 10, 2025
% Last Modified: June 10, 2025   

%   This code is an interpretation and adaptation of concepts from the from the
%   original MegAWES repository (Copyright 2021 Delft University of Technology),
%   which is also licensed under the Apache License, Version 2.0.


    filename = sprintf('%s_%d_%d.mp4', params_num.sim.aero_model, params_num.sim.path_type,params_num.sim.winding_number);
    % Set up the VideoWriter object
    animation_awe = VideoWriter(filename, 'Motion JPEG AVI');
    animation_awe.FrameRate = 15; % Adjust the frame rate as neede
    open(animation_awe);
    figure(99)
    set(figure(99),'units','normalized','outerposition',[0 0 1 1])

    for k = 1:length(t)
        % plot yx -> subplot(2,2,1)
        Visualisation(t(k), optimal_results(k).tether_pos(1,:)', optimal_results(k).tether_pos(2,:)', optimal_results(k).tether_pos(3,:)', optimal_results(k).pos_O, optimal_results(k).xn(9:11), params_num, h, 0, 90, 1)
        % plot xz -> subplot(2,2,2)
        Visualisation(t(k), optimal_results(k).tether_pos(1,:)', optimal_results(k).tether_pos(2,:)', optimal_results(k).tether_pos(3,:)', optimal_results(k).pos_O, optimal_results(k).xn(9:11), params_num, h, -90, 0, 2)
        % plot yz -> subplot(2,2,3)
        Visualisation(t(k), optimal_results(k).tether_pos(1,:)', optimal_results(k).tether_pos(2,:)', optimal_results(k).tether_pos(3,:)', optimal_results(k).pos_O, optimal_results(k).xn(9:11), params_num, h, 0, 0, 3)
        % plot xyz -> subplot(2,2,4)
        Visualisation(t(k), optimal_results(k).tether_pos(1,:)', optimal_results(k).tether_pos(2,:)', optimal_results(k).tether_pos(3,:)', optimal_results(k).pos_O, optimal_results(k).xn(9:11), params_num, h, -35, 15, 4)
            % Capture the current frame
        frame = getframe(gcf);
        % Write the frame to the video
        writeVideo(animation_awe, frame);
    end
    % Close the VideoWriter object
    close(animation_awe);
end