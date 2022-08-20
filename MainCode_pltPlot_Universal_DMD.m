%%pltPlot输出360动画通用程序，王佳琪，2020-12-21修改
%DMD论文再整理
%手动选择mat数据导入测试数据
%通用代码
%增加一判断失速功能，若失速，则自动进行双锁相；

clc
clear
close all
subfunction_path1='./subfunction/subfunction_1';
addpath(subfunction_path1);
[fname,location]=uigetfile({'*.mat';'*.*'},'mat参数文件读取','MultiSelect','on');%MultiSelect单选
load([location,'/','参数说明','/','parameter.mat']); %选择文件导入数据
disp(Note);
% % //======保存图像至指定文件夹===============  
save_directory = [strrep(location,'Database','DMD_DATA'),'叶顶流场动画结果输出-r-2-1',date];  %频谱图存储文件夹
if ~exist(save_directory)
    mkdir(save_directory)
else
    disp('文件夹存在！');
end
% % //========================================  
if isstr(fname)
   fname=cellstr(fname);
end
h=figure
c=colormap(jet(length(fname)+14));
sz = linspace(1,100,length(fname));
for i_file=1:length(fname)
Data = importdata(fullfile(location,char(fname(i_file)))); %选择文件导入数据
Data=V2Pa_Universal(Data,kulite_transform_ab);
Data(:,1:end-1)=Data(:,1:end-1);%-mean(Data(:,1:end-1));
[rpm,tsignal,xuhao,Rotor_Speed]=pltPlot_dB_universal(Data,fs,object,objectName,testTime,char(fname(i_file)),save_directory,0,'.mat');
% snapshots_compressor(rpm,tsignal,xuhao,save_directory,char(fname(i_file)));
%[tsignal2,EIGS]=computeKMD(rpm,tsignal,xuhao,save_directory,char(fname(i_file)));

[tsignal2,EIGS]=computeDMD(rpm,tsignal,xuhao,save_directory,char(fname(i_file)),sz(i_file),c(i_file,:));
%computeDMD_stream(rpm,tsignal,xuhao,save_directory,char(fname(i_file)),length(Data)/fs);
%[dmd_evals(:,i_file),tdmd_evals(:,i_file)]=computeDMD_streamtotal(rpm,tsignal,xuhao,save_directory,char(fname(i_file)),length(Data)/fs);

%[omega{i_file},EIGS{i_file},cond{i_file}]=computeDMD_test(Data,[4:16],fs,save_directory,char(fname(i_file)));
end
% figure;plot(real(cell2mat(omega)).');%衰减因子
% figure;plot(imag(cell2mat(omega)).');%频率--》相对速度
axis equal
axis([-1.1 1.1 -1.1 1.1]);
