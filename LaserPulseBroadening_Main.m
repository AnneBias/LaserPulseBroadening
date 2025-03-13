% Fig 4
% Fit of the FWHM of the simulated autocorrelation curve to CARPE measurements

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation and Fit of the FWHM of Simulated Autocorrelation to CARPE Measurements
%
% Description:
% This script simulates the broadening of a laser pulse due to material dispersion 
% and compares the simulated autocorrelation FWHM to experimentally measured data.
% The dispersion effects of the **microscope optics**, **fused silica windows**, 
% **ZnSe material**, and **biological tissue** are taken into account.
%
% Purpose:
% - Model the laser pulse propagation through dispersive media.
% - Compute the accumulated spectral phase due to different optical components.
% - Generate a simulated autocorrelation function.
% - Fit the simulated FWHM to experimental CARPE autocorrelation measurements.
%
% Methods:
% - Fourier transform-based pulse propagation.
% - Dispersion modeling using **Sellmeier equations** for refractive index calculations.
% - Numerical computation of spectral phase accumulation.
% - Second-order autocorrelation calculation.
% - Optimization of dispersion parameters (GDD, TOD) to match experiment.
%
% Input:
% - Measurement selection (e.g., 'Tibia' for biological tissue).
% - Experimental autocorrelation FWHM data.
% - Dispersion parameters (to be adjusted iteratively).
%
% Output:
% - Comparison between measured and simulated autocorrelation FWHM.
% - Visualization of the fit quality for different ZnSe thicknesses.
%
% Authors: [Scott Domingue, Anne Bias]
% Date: [Date]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
warning('off', 'MATLAB:polyfit:RepeatedPointsOrRescale')

%% --- Define Global Constants ---
global c                                % c is a global variable for convience
c =  299792458;                         % in m/s

% Time and frequency unit conversions
fs = 1e-15; 
ps = 1e-12;
nm = 1e-9;
THz = 1e12;

%% --- Define Fourier Axis Parameters --- 
NFT   =   2^12;                         % Number of points in the time-frequency domain
taxis = 5*ps;                           % Time axis span in seconds

% --- Time and Frequency Axes ---
tFT = (-NFT/2:NFT/2-1)/NFT*taxis;      
dt = tFT(2)-tFT(1);         
dnu = 1/(NFT*dt);                       
nu = (-NFT/2:NFT/2-1)*dnu;  
wFT =  2*pi*nu;                         % Zero-centered angular frequency axis (wFT = 0 == w0)

%% --- Define Initial Pulse Parameters ---
lambda0 = 1650 * 1e-9;                  % Central wavelength in meters
nu0 = c/lambda0;                        % Central frequency
w0 = nu0 * 2*pi;                        % Central angular frequency
lambda = c./(nu + nu0);                 % Wavelength array
lambdanm = lambda*1e9;                  
lambdanm(lambdanm<0) = nan;             % Remove negative values

% Define initial transform-limited pulse (Gaussian shape)
TauPulse = 65 * 1e-15;                  % Pulse duration (FWHM)
Et = exp(-2*log(2)*(tFT/TauPulse).^2);
ItTL = abs(Et).^2; ItTL = ItTL/max(ItTL);
Iw = abs( fftshift( fft( fftshift( Et ))) ).^2;
Iw = Iw/max(Iw);

%% --- Set Measurement Parameters ---
% Select measurement condition (tissue type)
Measurement = 'Tibia';
[MeasuredCarpeAC, MeasuredInsertion] = getMeasurementData(Measurement);

% Set dispersion parameters (to be optimized iteratively)
GDDlambda0Tissue = -859; 
GammaTODTissue = -3.7;
DeltaTau = 25;                          % Additional offset in fs


%% --- Calculate Microscope Dispersion (PhiBergamo) ---
x = 1400:50:1800;                       % Wavelength range for interpolation
GDDlambda0Microscope = -2013;
GammaTODMicroscope = -4.3;
LambdaDispersion0 = 1650;

% Interpolate microscope dispersion curve
y = GammaTODMicroscope*(x-LambdaDispersion0) + GDDlambda0Microscope;
GDDBergamo = interp1(x,y,lambdanm, [], 0); 
GDDBergamo(isnan(GDDBergamo)) = 0;

% Compute spectral phase of the microscope
PhiBergamo = GDDBergamo*1e-15^2/2.* wFT.^2 ; 

%% --- Calculate Tissue Dispersion (PhiTissue) ---
yTissue = GDDlambda0Tissue + GammaTODTissue*(x-LambdaDispersion0);
GDDTissue = interp1(x,yTissue,lambdanm, [], 0); 
GDDTissue(isnan(GDDTissue)) = 0;
PhiTissue = GDDTissue*1e-15^2/2.* wFT.^2 ; 

%% --- Compute Dispersion for Optical Components ---
% Material specifications
Material = 'ZnSe';   
FSLength = 6.34 * 1e-3;   % Fused silica window thickness in meters 

% Find valid wavelength indices to omit negative numbers in wavelength array
InIndx = findIndx(lambdanm, 2500):findIndx(lambdanm,500);
OutIndx = [1:InIndx(1)-1, InIndx(end)+1:NFT];

% Compute spectral phase for fused silica window
[PhiMat tmp] = matPhase(lambda,wFT,w0,InIndx,FSLength,'FS');
PhiMatExtra = PhiMat - PhiMat(NFT/2+1);

% Compute spectral phase for 1 mm of specified Material
[PhiMat GDD] = matPhase(lambda,wFT,w0,InIndx,1e-3,Material);
PhiMat = PhiMat - PhiMat(NFT/2+1);

% Compute material dispersion (GDD per mm)
p2 = polyfit(wFT(NFT/2+1+(-10:10)), PhiMat(NFT/2+1+(-10:10)), 2);
GDDMat = 2*p2(1)/fs^2;

%% --- Compute Total Dispersion and Autocorrelation ---
% Loop over all measured ZnSe thicknesses
for iL = 1:length(MeasuredInsertion)
    
    % Compute total spectral phase
    PhiC = PhiBergamo + PhiMat*MeasuredInsertion(iL) + PhiMatExtra + PhiTissue;

    % Compute intensity profile via inverse Fourier transform
    it = abs( ifftshift( ifft( ifftshift( sqrt(Iw) .* exp(-1i * PhiC ) ))) ).^2;
    it = it * trapz(ItTL)/trapz(it);
    
    % Compute electric field and FWHM duration
    Et = sqrt(it);
    FWHMDuration(iL) = fwhm(tFT/fs, it);
    
    % Compute numeric second harmonic autocorrelation
    indxRow = mod(reshape(0:NFT+1:(NFT+1)*NFT^2-1,NFT,NFT)+NFT^2/2,NFT^2)+1;
    Z = Et'*transpose(Et');                 % Outer product
    Z = fftshift(fft(Z(indxRow),[],1),1);   % Compute autocorrelation
    AC = sum(abs(Z).^2);
    AC = AC/max(AC); 
    
    % Compute autocorrelation FWHM with and without DeltaTau
    AutoCorrDuration(iL) = fwhm(tFT/fs, AC) + DeltaTau;       
    AutoCorrDuration_WithoutOffset(iL) = fwhm(tFT/fs, AC);    
end

%% --- Plot Results ---
figure(1);
clf;
hold on;

plot(MeasuredInsertion, AutoCorrDuration_WithoutOffset, 'k--', 'LineWidth', 1.5); % Simulated w/o DeltaTau
plot(MeasuredInsertion, AutoCorrDuration, 'k', 'LineWidth', 1.5); % Simulated w/ DeltaTau
plot(MeasuredInsertion, MeasuredCarpeAC, 'o-', 'LineWidth', 1.5, 'Color', [197/255, 55/255, 61/255]); % Measured data

legend({'Simulated w/o DeltaTau', 'Simulated w/ DeltaTau', 'Measured'});
xlabel('ZnSe Thickness [mm]', 'FontSize', 14);
ylabel('FWHM [fs]', 'FontSize', 14);
ylim([90 170]);
xlim([3 12]);

return