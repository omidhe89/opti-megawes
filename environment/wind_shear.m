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

function MeanWind_W = wind_shear(DCM, altitude, params_sym)
    % WIND_SHEAR Calculates the mean wind vector in the World frame, considering wind shear.
    %
    % This function computes the local mean wind velocity at a given altitude,
    % accounting for a logarithmic wind shear profile. It then transforms this
    % velocity from a local wind-aligned frame to the global World frame.
    %
    % Inputs:
    %   DCM         - Rotation matrix (3x3) from the World frame to the body reference frame.
    %                
    %   altitude    - Scalar, The altitude (height above ground) at which to calculate
    %                 the wind speed. Please note that in some conventions, altitude
    %                 might be `-Pos_O(3)` if `Pos_O(3)` is the Z-coordinate in a
    %                 North-East-Down (NED) or similar frame. Ensure `altitude` is positive upwards.
    %   params_sym  - struct, Symbolic parameters.
    %
    % Outputs:
    %   MeanWind_W  - Vector (3x1), The mean wind velocity vector at the specified
    %                 altitude, expressed in the World (global) coordinate frame [vx_W; vy_W; vz_W].

    % Extract base wind speed and direction from symbolic parameters.
    base_wind_speed = params_sym.env.wind.base_velocity;
    wind_direction = params_sym.env.wind.direction;

   
    BaseWind = -base_wind_speed .* [cos(wind_direction); sin(wind_direction); 0];

    % Apply the logarithmic wind shear model to get the wind speed at the given altitude.
    % The `log_wind` function should compute the wind velocity magnitude based on altitude
    % and then potentially orient it.
    % Note: The `log_wind` function signature usually takes `altitude` and `base_wind_speed`
    Windspeed_altitude = log_wind(altitude, BaseWind, params_sym);
    MeanWind_O = -1 .* (DCM' * Windspeed_altitude);

    % Apply a final transformation from the intermediate 'O' frame to the World 'W' frame
    MeanWind_W = transformFromOtoW(wind_direction, MeanWind_O);
end