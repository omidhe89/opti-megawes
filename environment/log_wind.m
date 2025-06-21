% Copyright (c) 2025 [Your Name or Organization Name]
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

function v_w_p = log_wind(z_pos, base_wind_speed, params_sym)
    % LOG_WIND Calculates wind speed at a given height using a logarithmic wind profile.
    %
    % This function implements the logarithmic wind shear law, which is commonly
    % used to model how wind speed varies with altitude in the atmospheric
    % boundary layer over homogeneous terrain.
    %
    % Inputs:
    %   z_pos           - Scalar, The height above the ground (altitude) at which
    %                     to calculate the wind speed (m). This corresponds to `z` in the formula.
    %   base_wind_speed - Scalar, The known wind speed at a reference height (m/s).
    %                     This corresponds to `v_ref` in the formula.
    %   params_sym      - struct, Symbolic parameters.
    %
    % Outputs:
    %   v_w_p           - Scalar, The calculated wind speed at the specified height `z_pos` (m/s).
    %                     This corresponds to `v(z)` in the formula.
    %

    z0    = params_sym.env.wind.z0;    % Surface roughness length (m)
    z_ref = params_sym.env.wind.z_ref; % Reference height (m)

    % Calculate the wind speed at the given height (z_pos) using the logarithmic profile.
    % Make sure z_pos is not zero or less than z0 to avoid log(0) or log(negative) issues.
    v_w_p = base_wind_speed .* (log10(z_pos / z0) / log10(z_ref / z0));

end