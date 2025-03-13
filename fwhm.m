function [fwhm, i0] = fwhm(t, It)
% FWHM calculates the Full-Width-Half-Maximum (FWHM) of a given intensity profile.
%
% This function normalizes the input intensity profile and finds the 50% 
% intensity points on either side of the peak to determine the FWHM. 
% The function does not subtract any background (noise).   
% 
% Arguments:
%   t       - Vector containing time values
%   It      - Intensity profile (does not need to be pre-normalized)
%
% Returns:
%   fwhm    - Full-Width-Half-Maximum of the intensity profile
%   i0      - Index of the center of the FWHM range

%% --- Parameter Initialization & Normalization ---
allowed_err = 0.001;                    % Tolerance for detecting the 50% threshold
It = It ./ max(It);                     % Normalize intensity to its max value

%% --- Locate the Left Half-Maximum Point ---
[temp, imax] = max(It);                 % Find the peak index
t1 = 1;
t2 = imax;

% Search for the left 50% intensity point
counter = 1;
while counter < imax
    t3 = round((t2+t1)/2);
    err = (It(t3)-0.5)/0.5;
    if t2 == t1 + 1
        break 
    elseif It(t3)>0.5
        t2 = t3;
    else 
        t1 = t3;
    end
    counter = counter+1;
end
i_left = t1; 


% Linear interpolation for a more precise half-maximum point
if (0.5 - It(t1)) < allowed_err
    t_left = t(t1);
else
    n = round((It(t2) - It(t1))/0.5/allowed_err);
    t_ = linspace(t(t1), t(t2), n);
    It_ = interp1(t(t1:t2), It(t1:t2), t_);
    
    err = 1; t1 = 1;  t2 = n;
    while abs(err)>allowed_err+0.001
        t3 = round((t2+t1)/2);
        err = (It_(t3)-0.5)/0.5;
        if err > 0
            t2 = t3;
        else
            t1 = t3;
        end
    end
    t_left = t_(t3);
end

%% --- Locate the Right Half-Maximum Point ---
t1 = imax; 
t2 = length(t);

% Search for the right 50% intensity point
counter = 1;
while counter < length(t)-imax
    t3 = round((t2+t1)/2);
    err = abs(It(t3)-0.5)/0.5;
    if t2 == t1 + 1
        break 
    elseif It(t3)>0.5
        t1 = t3;
    else 
        t2 = t3;
    end
    counter = counter + 1;
end
i_right = t1;

% Linear interpolation for a more precise half-maximum point
if (0.5-It(t1)) < allowed_err
    t_right = t(t1);
else
    n = round((It(t1) - It(t2))/0.5/allowed_err);
    t_ = linspace(t(t1), t(t2), n);
    It_ = interp1(t(t1:t2), It(t1:t2), t_);
    
    err = 1; t1 = 1;  t2 = n;
    while abs(err)>allowed_err+0.001
        t3 = round((t2+t1)/2);
        err = (It_(t3)-0.5)/0.5;
        if err > 0
            t1 = t3;
        else
            t2 = t3;
        end
    end
    t_right = t_(t3);
end

%% --- Final Calculation of FWHM and Center Index ---
fwhm = abs(t_right - t_left);
i0 = round((i_right + i_left)/2);