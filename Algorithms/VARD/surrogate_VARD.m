% Code written by Y. Kaganovsky,  Dec. 2014
% This function computes the 1D surrogate for the variance S_{v_j} used by VARD
% as well as its 1st and 2nd derivatives

function [fun der1 der2] = surrogate(v,v_prev,sbp_E,C,Z)
fun = sbp_E.*exp(Z.*(v-v_prev))./Z./2-log(v)./2+C.*v./2;
der1 = sbp_E.*exp(Z.*(v-v_prev))./2-1./v./2+C./2;
der2 = sbp_E.*Z.*exp(Z.*(v-v_prev))./2+1./(v.^2)./2;
end