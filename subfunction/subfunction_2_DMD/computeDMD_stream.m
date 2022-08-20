%请参考https://github.com/cwrowley/dmdtools
function [omega,EIGS,cond]=computeDMD_stream(rpm,tsignal,xuhao,save_directory,name,T)
dt=T/length(rpm);

for kk=1:length(rpm)
    len(kk)=length(rpm{1,kk});
end
[q1,q2]=min(len)
sum=zeros(min(len),10);
for kk=1:length(rpm)
    sum=sum+rpm{1,kk}(1:min(len),:);
end
av=sum/length(rpm);

for kk=1:length(rpm)
    VORTALL(:,kk)=reshape(rpm{1,kk}(1:min(len),:)-av,1,[])';%按列
end

%%
% set example parameters here
max_rank = 10;      % maximum allowable rank of the DMD operator
                    %   (set to zero for unlimited)
m = size(VORTALL,2)-1;            % total number of snapshots to be processed
n = size(VORTALL,1);           % number of states
%noise_cov = 1e-4;   % measurement noise covariance
%rng(0)              % seed the random number generator
streaming = true;   % true=use streaming DMD, false=use batch-processed DMD
%%


% X = VORTALL(:,1:end-1);
% X2 = VORTALL(:,2:end);
% [U,S,V] = svd(X,'econ');
% cond=sqrt(max(diag(S))/min(diag(S))) %条件数越大，矩阵病态越验证

% %%  Compute DMD (Phi are eigenvectors)
% r =30;%length(rpm)-1;  % truncate at 21 modes
% U = U(:,1:r);
% S = S(1:r,1:r);
% V = V(:,1:r);
% Atilde = U'*X2*V*inv(S);
% [W,eigs] = eig(Atilde);
% Phi =real(X2*V*inv(S)*W);
% EIGS=diag(eigs);
% omega = log(EIGS)/dt; % continuous-time eigenvalues
% %abs(omega)



%% collect the data

disp('Collecting data')
tic
if ~streaming
    % standard algorithm: batch-processing
    X = zeros(n, m);
    Y = zeros(n, m);
    yk = VORTALL(:,1);
    for k = 1:m
        xk = yk;
        yk = VORTALL(:,k+1);
        X(:,k) = xk;
        Y(:,k) = yk;
    end
    [Qx, S, V] = svd(X, 0);
    Ktilde = Qx' * Y * V * pinv(S);
else
    % streaming algorithm
    sdmd = StreamingDMD(max_rank);
    yk = VORTALL(:,1);
    for k = 1:m
        xk = yk;
        yk = VORTALL(:,k+1);
        sdmd = sdmd.update(xk, yk);
    end
end
elapsed_time = toc;
fprintf('  Elapsed time: %f seconds\n', elapsed_time)




%% compute DMD spectrum
disp('Computing spectrum')
tic
if ~streaming
    % standard DMD spectrum
    [evecK, evals] = eig(Ktilde);
    evals = diag(evals);
    modes = Qx * evecK;
else
    % streaming algorithm
    [modes, evals] = sdmd.compute_modes();
end
elapsed_time = toc;
fprintf('  Elapsed time: %f seconds\n', elapsed_time)

% calculate corresponding frequencies
fdmd = abs(angle(evals)) ./ (2 * pi * dt);
ydmd = zeros(length(fdmd),1);
for ii = 1:length(fdmd)
    ydmd(ii) = norm(modes(:,ii)).*abs(evals(ii));
end
ydmd = ydmd./max(ydmd);

figure(1)
stem(fdmd,ydmd,'o-')
xlabel('Frequency')
ylabel('Magnitude')








% 
% 
% 
% tsignal2.Nvar= tsignal.Nvar;
% tsignal2.varnames= tsignal.varnames;
% for kk=1:r
% tsignal2.surfaces(kk).zonename= tsignal.surfaces(q2(1)).zonename;
% tsignal2.surfaces(kk).x= tsignal.surfaces(q2(1)).x;
% tsignal2.surfaces(kk).y= tsignal.surfaces(q2(1)).y;
% tsignal2.surfaces(kk).z= tsignal.surfaces(q2(1)).z;
% tsignal2.surfaces(kk).v= reshape(Phi(:,kk),1,min(len),10); 
% %tsignal2.surfaces(kk).v(1,xuhao,:)=1.1*max(Phi(:,kk));
% tsignal2.surfaces(kk).solutiontime=kk;
% end
% mat2tecplot(tsignal2,[save_directory,'\',name,'.plt']);
% %% Plot DMD modes
% 
% % for i=10:2:20
% %     figure
% %     imagesc(reshape(real(Phi(:,i)),min(len),10)); % plot vorticity field
% %     figure
% %     imagesc(reshape(real(Phi(:,i)),min(len),10)); % plot vorticity field
% % % 
% % %     plotCylinder(reshape(real(Phi(:,i)),min(len),10),min(len),10);
% % %     plotCylinder(reshape(imag(Phi(:,i)),min(len),10),min(len),10);
% % end
% % 
% %%  Plot DMD spectrum
% h=figure
% theta = (0:1:100)*2*pi/100;
% plot(cos(theta),sin(theta),'k--') % plot unit circle
% hold on, grid on
% real_eigs=real(diag(eigs));
% imag_eigs=imag(diag(eigs));
% scatter(real_eigs,imag_eigs,'ok')
% for i=1:length(eigs)
% text(real_eigs(i),imag_eigs(i),num2str(i))
% end
% axis([-1.1 1.1 -1.1 1.1]);
saveas(h,[save_directory,'\',name,'.png'])
end