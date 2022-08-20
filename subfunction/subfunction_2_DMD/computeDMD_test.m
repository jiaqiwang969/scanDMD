
function [omega,EIGS,cond]=computeDMD_test(Data,object,fs,save_directory,name)
 dt=1/fs;
% 
% for kk=1:length(rpm)
%     len(kk)=length(rpm{1,kk});
% end
% [q1,q2]=min(len)
% for kk=1:length(rpm)
%     VORTALL(:,kk)=reshape(rpm{1,kk}(1:min(len),:),1,[])';%按列
% end

X = Data(object,1:end-1);
X2 = Data(object,2:end);
[U,S,V] = svd(X,'econ');
cond=sqrt(max(diag(S))/min(diag(S))) %条件数越大，矩阵病态越验证

%%  Compute DMD (Phi are eigenvectors)
r =11; %length(rpm)-1;  % truncate at 21 modes
U = U(:,1:r);
S = S(1:r,1:r);
V = V(:,1:r);
Atilde = U'*X2*V*inv(S);
[W,eigs] = eig(Atilde);
Phi =real(X2*V*inv(S)*W);
EIGS=diag(eigs);
omega = log(EIGS)/dt; % continuous-time eigenvalues
%abs(omega)





%% Plot DMD modes

% for i=10:2:20
%     figure
%     imagesc(reshape(real(Phi(:,i)),min(len),10)); % plot vorticity field
%     figure
%     imagesc(reshape(real(Phi(:,i)),min(len),10)); % plot vorticity field
% % 
% %     plotCylinder(reshape(real(Phi(:,i)),min(len),10),min(len),10);
% %     plotCylinder(reshape(imag(Phi(:,i)),min(len),10),min(len),10);
% end
% 
%%  Plot DMD spectrum
h=figure
theta = (0:1:100)*2*pi/100;
plot(cos(theta),sin(theta),'k--') % plot unit circle
hold on, grid on
real_eigs=real(diag(eigs));
imag_eigs=imag(diag(eigs));
scatter(real_eigs,imag_eigs,'ok')
for i=1:length(eigs)
text(real_eigs(i),imag_eigs(i),num2str(i))
end
axis([-1.1 1.1 -1.1 1.1]);
saveas(h,[save_directory,'\',name,'.png'])
end