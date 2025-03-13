function [MeasuredCarpeAC, MeasuredInsertion] = getMeasurementData(Measurement)
% GETMEASUREMENTDATA retrieves experimental FWHM measurement data for different tissue types.
%
% This function returns the measured FWHMs (in fs) of the second-order 
% interferometric autocorrelation of the laser pulse (measured with CARPE, APE) 
% as well as the corresponding ZnSe thicknesses (in mm) used for dispersion correction.
%
% Arguments:
%   Measurement       - String specifying the measurement condition ('Microscope', 'Tibia', etc.)
%
% Returns:
%   MeasuredCarpeAC   - Measured FWHM values (fs) of the autocorrelation measurement
%   MeasuredInsertion - ZnSe thicknesses (mm) for dispersion correction

    switch Measurement

        case 'Microscope'                      
            % Measurement taken in air under the objective for the 
            % iterative determination of microscope dispersion parameters 
            % GDDlambda0Microscope and GammaTODMicroscope.
            %
            % --- Important: This simulation must be performed first!
            % The iteratively determined microscope parameters
            % (GDDlambda0Microscope, GammaTODMicroscope) serve as the basis
            % for all subsequent tissue simulations.
            %    
            % --- Setup Instructions:
            % - set PhiMicroscope = 0 in the main script 
            % - deactivate Kerr pulse broadening & convolution (itPulse == it)
            % -> The actual PhiMicroscope calculation will be performed in PhiTissue 
            %    (GDDlambda0Tissue, select GammaTODTissue, DeltaTau = 0).
            %
            % Once the microscope parameters are determined iteratively: 
            % - insert them in the "Calculation of Microscope Dispersion 
            %   (PhiBergamo)" section in GDDlambda0Microscope and GammaTODMicroscope
            % - reactivate Kerr pulse broadening & convolution
            % - continue with tissue simulations to determine their 
            %   GDDlambda0Tissue, select GammaTODTissue, DeltaTau

            MeasuredCarpeAC = [117 111 104 99 96 94 92 92 93 96 99 103 107 112 118 124 132 140 145]; 
            MeasuredInsertion = (3:0.5:12);    

        case 'Tibia'                       
            MeasuredCarpeAC = [159 154 145 140 134 129 125 123 122 120 120 121 124 127 130 132 137]; 
            MeasuredInsertion = (4:0.5:12);    

        otherwise
            error('Unknown measurement: %s.', Measurement);

    end
end