function n = sellmeier(lambda,material)
% SELLMEIER computes the refractive index for a given material and 
% wavelength using the Sellmeier equation.
%
% Arguments:
%   lambda   - Wavelength(s) in microns (must be provided in µm)
%   material - String specifying the material 
%
% Returns:
%   n        - Refractive index for the given material at the specified wavelength(s)
%
% Available materials:
%   - 'FS'   (Fused Silica)  - added 08/23/2010    ref 1
%   - 'ZnSe' (Zinc Selenide) - added 08/24/2010    ref 2
%
% References:
%   ref 1) http://www.cvimellesgriot.com/products/Documents/Catalog/Dispersion_Equations.pdf
%   ref 2) tatian_AO_1984
%
% Version I - sdomingue - 08/23/2010

%% --- Initialize Sellmeier Coefficients ---
b = zeros(1,3); 
c = zeros(1,3);
n = zeros(1, length(lambda));

%% --- Assign Coefficients Based on Material ---
switch material
    case 'FS'
        N = 1;
        b = [0.6961663, 0.4079426, 0.8974794];
        c = [0.00467914826, 0.0135120631, 97.9340025];

    case 'ZnSe'
        N = 2;
        b = [4.45813734, 0.467216334, 2.8956629];
        c = [0.200859853, 0.391371166, 47.1362108];

    otherwise
        error('Unknown material: %s.', material);
            
end

%% --- Compute Refractive Index using Sellmeier Equation ---
n = sqrt(1 + ...
    b(1)*lambda.^2./(lambda.^2 - c(1)^N) + ...
    b(2)*lambda.^2./(lambda.^2 - c(2)^N) +  ...
    b(3)*lambda.^2./(lambda.^2 - c(3)^N));