function [ output_curve ] = dialogdiscard( message,curve )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
choice = questdlg(message, ...
	'Please take note of this decision', ...
	'Yes','NO, I will discard this cell','NO, I will discard this cell');
% Handle response
switch choice
    case 'Yes'
        disp('this cell will be kept')
        output_curve=curve;
    case 'NO, I will discard this cell'
        disp([' This is wiser'])
       output_curve=[];
   
end

end

