function [a] = fun_gastransfer_a(par_gastransfer_a,a_default)
% fun_gastransfer_a
%
%   ***********************************************************************
%   *** air-sea gas transfer coefficient for the config file **************
%   ***********************************************************************
%
%   par_gastransfer_a, if set (not empty) in the parameter file, replaces
%   the coefficient calibrated for the input format and wind-speed product
%   (a_default)
%
%   26/10/16: created
%
if isempty(par_gastransfer_a)
    a = a_default;
else
    a = par_gastransfer_a;
end
