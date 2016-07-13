% This file includes a demo of the VARD algorithm discussed in the paper 
%``Alternating Minimization Algorithm with Automatic Relavence Determination for Transmission Tomography under Poisson Noise
% by Y. Kaganovsky, S. Han, S. Degirmenci, D. Politte, D. Brady, J. O'Sullivan. and L. Carin.
% http://arxiv.org/pdf/1412.8464v1.pdf

% Version 1.0, Dec. 2014
% Code written by Y. Kaganovsky  and S. Degirmenci

% In order to make the code suitable for any personal computer we
% have decreased the image resolution from 256x256 to 128x128
% the rest of the parameters and initializations are the same as in the paper


clear all
close all
clc

addpath('Algorithms/VARD');
addpath('System Matrix');

fprintf('This script will perform reconstruction using only the VARD algorithm \n');


% Parameters %
mu_water =0.0190;                  % linear attenuation coeffecient of water in mm^{-1}
eta=1e5;                           % Poisson rate for air scan
img_size=128;                      % number of pixels along each dimension of the image

% plot truth %
object=phantom(img_size);          % Ground truth for attenuation image
fig=figure;
imagesc(object);
colormap gray
colorbar;
title('Ground Truth', 'fontsize', 14);
caxis([0.15 0.35]);
set(gca, 'fontsize', 14) ;
set(fig, 'units', 'inches', 'position', [4 1 8 7]);
fprintf('You can examine the figure and then press any key to continue to reconstruction \n');
pause

%%%%%%%%%%%%%  Load Data and Matrices to Memory %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Loading data to memory (Initialization step) \n');

%% load system matrix %
fprintf('Loading system matrix... ');
load(['H_128pix_full_offset_synth.mat'],'H');                     % Here the system matrix is H (in the manuscript we denote it by Phi)
fprintf('Done \n')

%% Generate data using Poisson random realization   %%
fprintf('Generating random realizations for Poisson data... ');
Proj=H.*mu_water*object(:);        	                              % Projections
y=poissrnd(eta.*exp(-Proj));                                      % Poisson Data 
fprintf('Done \n')

% Find circular domain of interest  %
[rr cc] = meshgrid(1:img_size);
center=round(img_size./2);
ind_circ = find(sqrt((rr-center).^2+(cc-center).^2)<=center);
N=numel(ind_circ);                                                % Number of pixels in the circular domain


% Keep parts of system matrix and measurements that are relavent to the circular domain of interest %
H=H(:,ind_circ).*mu_water;
ind_nonempty=find(sum(H,2)~=0);
H=H(ind_nonempty,:);
y=y(ind_nonempty);

% Calculate transpose of system matrix %
fprintf('Calculating transpose of system matrix... ');
H_t=H.';
fprintf('Done \n');

% Calculate squared system matrix %
fprintf('Calculating squared system matrix... ');
H_sq=H.^2;
fprintf('Done \n')

% Calculate transpose of squared system matrix %
fprintf('Calculating transpose of squared system matrix... ');
H_sq_t=H_sq.';
fprintf('Done \n')

%%%%%% create prior precision matrix %%%%%%%%%
fprintf('Generating neighborhood matrix... ');
Psi_h=speye(img_size.^2)+sparse(diag(-ones(img_size.^2-1,1),1));                    % allocate memory for a sparse matrix
Psi_v=speye(img_size.^2)+sparse(diag(-ones(img_size.^2-img_size,1),-img_size));       
Psi_h=Psi_h(ind_circ,ind_circ);                                                     % difference operator for horizontal neighbors  
Psi_v=Psi_v(ind_circ,ind_circ);                                                     % difference operator for vertical neighbors  
Psi=[Psi_h; Psi_v];                                           
fprintf('Done \n')

fprintf('Calculating backprojected data... ');

BP_y=H.'*y;
Z_ml=max(sum(H,2));                        % For MLE and MAP
Z_vard=max(sum(H+0.5.*H.^2,2));            % For VARD
fprintf('Done \n')

fprintf('Initialization finished \n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N_iter=2000;                                                          % number of iterations
% reconstruct using VARD   (NO TUNING PARAMETERS NEEDED)
fprintf('Starting Variational Automatic Relevance Determination (VARD) \n');
VARD;      %  we use the overcomplete sparse representation given by Eq. (6.1) in the manuscript

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('The demo is finished \n');
