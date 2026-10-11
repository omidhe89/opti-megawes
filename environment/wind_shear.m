% Copyright (c) 2025 Omid Heydarnia
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%
%    http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
%
% Author: Omid Heydarnia
% Email: omid.heydarnia@ugent.be
% Date: June 12, 2025
% Last Modified: June 12, 2025

function MeanWind_B = wind_shear(DCM, pos_W, params_sym)
    % WIND_SHEAR Calculates the mean wind vector in the World frame, considering wind shear.
    %
    % This function computes the local mean wind velocity at a given altitude,
    % accounting for a logarithmic wind shear profile. It then transforms this
    % velocity from a local wind-aligned frame to the global World frame.
    %
    % Inputs:
    %   DCM         - Rotation matrix (3x3) from the World frame to the body reference frame.
    %                
    %   pos_W       - Vector (3x1), Position in the World frame. The Z-coordinate (pos_W(3))
    %                 represents the altitude (height above ground).
    %   params_sym  - struct, Symbolic parameters containing wind properties (base velocity,
    %                 direction, shear parameters, etc.).

    %
    % Outputs:
    %   MeanWind_B  - Vector (3x1), The mean wind velocity vector at the specified
    %                 altitude, expressed in the Body (local) coordinate frame [vx_B; vy_B; vz_B].

    % Extract base wind speed and direction from symbolic parameters.
    base_wind_speed = params_sym.env.wind.base_velocity;
    wind_direction = params_sym.env.wind.direction;
    altitude = pos_W(3); % Altitude of the kite (assuming pos_W(3) = -pos_O(3) is Z-coordinate in ground frame)


    % Apply the logarithmic wind shear model to get the wind speed at the given altitude.
    Windspeed_altitude = log_wind(altitude, base_wind_speed, params_sym);
    Windspeed = Windspeed_altitude .* [cos(wind_direction); sin(wind_direction); 0];
    MeanWind_B = DCM' * Windspeed;
end