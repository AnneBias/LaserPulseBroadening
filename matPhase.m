function [Phi GDD TOD] = matPhase(lambda, wFT, w0, InIndx, L, MAT)
% MATPHASE calculates the phase accumulation of an optical field 
% through a material using the Sellmeier equation and polynomial fitting.
%
% Arguments:
%   lambda   - Wavelengths 
%   wFT      - Frequency grid (angular frequency deviations from w0)
%   w0       - Central angular frequency
%   InIndx   - Indices of the valid wavelength range for material calculation
%   L        - Material thickness
%   MAT      - String specifying the material (e.g., 'FS', 'ZnSe')
%
% Returns:
%   Phi      - Computed spectral phase (relative to the group velocity)
%   GDD      - Group delay dispersion 
%   TOD      - Third-order dispersion 
%
% Global Variables:
%   c        - Speed of light (must be defined in the main script)
%
% Notes:
% - This function determines the refractive index using the Sellmeier equation.
% - The spectral phase is computed by subtracting the group velocity component.
% - GDD and TOD are extracted using polynomial fitting.

global c;

%% --- Initialization ---
NFT = length(wFT);
OutIndx = [1:(InIndx(1)-1), (InIndx(end)+1):length(lambda)];

%% --- Compute Refractive Index using Sellmeier Equation ---
% Compute refractive index only for the valid wavelength range
n(InIndx) = sellmeier(lambda(InIndx)./1e-6, MAT);
% For wavelengths outside the valid range, use the mean of the known values
n(OutIndx) = mean(n(InIndx)); 
% Note: This could alternatively be set to zero, but caution is needed 
% when integrating this into a broader spectrum simulation.

%% --- Compute Wave Vector and Spectral Phase ---
% Wave vector representation (β = n * ω / c) is the most general expression 
% for the phase accumulation of an optical field through some material
Beta =  n.*(wFT+w0)/c ;

% Polynomial fit range (centered around w0)
LengthFit = 100;   
% Restrict the polynomial fitting range to improve accuracy at the central frequency.
% Increasing this range averages over nearby frequencies, affecting accuracy.
    
%% --- Compute Group Velocity and Remove Linear Component ---
GV = polyfit(wFT(round(NFT/2)+1+(-LengthFit:LengthFit)), ...
    Beta(round(NFT/2)+1+(-LengthFit:LengthFit)),1);
GV = polyval(GV, wFT);  
% Fit a polynomial over the restricted range but extrapolate linearly over 
% the full frequency domain.
    
Beta2 = Beta - GV;          % Remove linear phase contribution (group velocity)
Beta2(OutIndx) = 0;         % Ensure values outside range are set to zero

%% --- Compute Final Spectral Phase ---
Phi  = Beta2*L;             % Apply material thickness
Phi = Phi-min(Phi(InIndx)); % Remove phase offset
    
%% --- Extract Dispersion Parameters (GDD & TOD) ---
% GDD (fs²) - Second derivative of spectral phase
GDD = polyfit(wFT(round(NFT/2)+1+(-LengthFit:LengthFit)),...
    Phi(round(NFT/2)+1+(-LengthFit:LengthFit)),2);
GDD = GDD(1) * 2 / 1e-15^2;
    
% TOD (fs³) - Third derivative of spectral phase
TOD = polyfit(wFT(round(NFT/2)+1+(-LengthFit:LengthFit)),...
    Phi(round(NFT/2)+1+(-LengthFit:LengthFit)), 3);
TOD = TOD(1) * 6 / 1e-15^3;

return