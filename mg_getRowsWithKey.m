% Mex function, everything is commented
%  Author: Matthew Gunn
%  Updated:    11/10/15
%
%  Usage:      [key_unique, keymap] = mg_getRowsWithKey(keyvals)
%
%             This function 
%
%  keyvals:   This can either be (i)  a vector of real doubles
%                             or (ii) a cell array of strings
%
%  key_unique: This contains only the UNIQUE members of keyvals
%
%  keymap:     A cell array where keymap{i} will hold a vector containing
%              all the rows of keyvals that are equivalent to key_unique(i)
%
%  Example usage:
%    [key_unique, keymap] = mg_getRowsWithKey(group_id_vector)
% 
%    for(i=1:length(key_unique)),
%       group_i_key  = key_unique(i);
%       group_i_data = data(keymap{i},:);
%    end
