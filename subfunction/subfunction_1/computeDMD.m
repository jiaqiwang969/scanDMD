
function [tsignal2,EIGS]=computeDMD(rpm,tsignal,xuhao,save_directory,name,sz,c)
% dt=T/length(rpm);

for kk=1:length(rpm)
    len(kk)=length(rpm{1,kk});
end
[q1,q2]=min(len)
for kk=1:length(rpm)
    VORTALL(:,kk)=reshape(rpm{1,kk}(1:min(len),:),1,[])';%°´ÁÐ
end

X = VORTALL(:,1:end-1);
X2 = VORTALL(:,2:end);
[U,S,V] = svd(X,'econ');

%%  Compute DMD (Phi are eigenvectors)
r =11;  % truncate at 21 modes
U = U(:,1:r);
S = S(1:r,1:r);
V = V(:,1:r);
Atilde = U'*X2*V*inv(S);
[W,eigs] = eig(Atilde);
Phi =real(X2*V*inv(S)*W);
EIGS=diag(eigs);



tsignal2.Nvar= tsignal.Nvar;
tsignal2.varnames= tsignal.varnames;
for kk=1:r
tsignal2.surfaces(kk).zonename= tsignal.surfaces(q2(1)).zonename;
tsignal2.surfaces(kk).x= tsignal.surfaces(q2(1)).x;
tsignal2.surfaces(kk).y= tsignal.surfaces(q2(1)).y;
tsignal2.surfaces(kk).z= tsignal.surfaces(q2(1)).z;
tsignal2.surfaces(kk).v= reshape(Phi(:,kk),1,min(len),10); 
%tsignal2.surfaces(kk).v(1,xuhao,:)=1.1*max(Phi(:,kk));
tsignal2.surfaces(kk).solutiontime=kk;
end
mat2tecplot(tsignal2,[save_directory,'\',name,'.plt']);
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

theta = (0:1:100)*2*pi/100;
plot(cos(theta),sin(theta),'k--') % plot unit circle
plot(0.5*cos(theta),0.5*sin(theta),'k--') % plot unit circle
hold on, grid on
real_eigs=real(diag(eigs));
imag_eigs=imag(diag(eigs));
scatter(real_eigs,imag_eigs,sz,'MarkerFaceColor',c,'MarkerEdgeColor',c,'MarkerFaceAlpha',.8,'MarkerEdgeAlpha',.8)
hold on

% for i=1:length(eigs)
% text(real_eigs(i),imag_eigs(i),num2str(i))
% end

%saveas(h,[save_directory,'\',name,'.png'])
end
