function  indx = findIndx (input1, input2)
% FINDINDX finds the index of the closest value in input1 to input2.
%
% Arguments:
%   input1 - Vector of values
%   input2 - Reference value to compare against
%
% Returns:
%   indx   - Index of the closest value in input1 to input2
%
% Copyright (C) 2016 sdomingue
% Author: sdomingue <sdomingue@KML55>
% Created: 2016-08-11

    [tmp TMP] = sort(abs( input1 - input2 ));
    indx = TMP(1);
return
