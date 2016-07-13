% The algorithm implemented here belongs to the class of algorithms presented in the paper
% ``Alternating Minimization Algorithm for Transmission Tomography'', by O'Sullivan and Benac, IEEE Trans. on Medical Imaging, 2007
% Code written by Y. Kaganovsky,  Dec. 2014

x = zeros(N,1);                       % initial guess for posterior mean
obj_fun = zeros(1,N_iter);            % vector to save the history of the objective

reverseStr  =  '';
tic

for iter = 1:N_iter
    msg  =  sprintf('Now running iteration No. %d/%d', iter,N_iter);
    fprintf([reverseStr, msg]);
    reverseStr  =  repmat(sprintf('\b'), 1, length(msg));
    proj = H_t.'*x;                                                     % (line 4 of Algorithm 3 in the manuscript),  this is faster than H*x due to the way Matlab stores the matrix in memory
    E = eta.*exp(-proj);                                                % (line 5 of Algorithm 3 in the manuscript)
    BP_E = H.'*E;                                                       % (line 6 of Algorithm 3 in the manuscript)
    %Note that Matlab automatically exploits multiple CPU cores and parallel computing when performing the vectorized operations below
    x = x+log(BP_E./BP_y)./Z_ml;                                        % (line 8 of Algorithm 3 in the manuscript)
    x(find(x<0)) = 0;                                                   % Threshold negative values
    obj_fun(iter) = y.'*proj+sum(E);                                    % evaluate objective at current iteration
end
time = toc;
fprintf(['\nFinished! total run time was ' num2str(round(time./60)) ' mins \n']);      % display total run time

%%%  plot figures     %%%
% plot reconstructed image %
x_plot = zeros(img_size);
x_plot(ind_circ) = x;
fig = figure;
imagesc(x_plot);
colormap gray
caxis([0.15 0.35]);
colorbar;
title(['Reconstructed Image | MLE after '  num2str(N_iter) ' Iterations'],'fontsize', 14);
set(gca, 'fontsize', 14)
set(fig, 'units', 'inches', 'position', [4 1 8 7]);


% plot objective history %
fig = figure;
semilogx(1:iter,obj_fun,'linewidth',3);
title('Objective vs Iteration Number | MLE','fontsize', 14);
set(gca, 'fontsize', 14)
xlabel('Iteration Number (log-scale)');
ylabel('Objective function value');
xlim([0 N_iter]);
grid on
set(fig, 'units', 'inches', 'position', [4 1 8 7]);


RMSE_vard = norm(x-object(ind_circ))./norm(object).*100;
fprintf(['The RMSE for MLE is ' num2str(RMSE_vard) '%% \n']);
