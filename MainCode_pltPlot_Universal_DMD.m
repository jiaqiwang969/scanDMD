%%pltPlot���360����ͨ�ó�����������2020-12-21�޸�
%DMD����������
%�ֶ�ѡ��mat���ݵ����������
%ͨ�ô���
%����һ�ж�ʧ�ٹ��ܣ���ʧ�٣����Զ�����˫���ࣻ

clc
clear
close all
subfunction_path1='./subfunction/subfunction_1';
addpath(subfunction_path1);
[fname,location]=uigetfile({'*.mat';'*.*'},'mat�����ļ���ȡ','MultiSelect','on');%MultiSelect��ѡ
load([location,'/','����˵��','/','parameter.mat']); %ѡ���ļ���������
disp(Note);
% % //======����ͼ����ָ���ļ���===============  
save_directory = [strrep(location,'Database','DMD_DATA'),'Ҷ����������������-r-2-1',date];  %Ƶ��ͼ�洢�ļ���
if ~exist(save_directory)
    mkdir(save_directory)
else
    disp('�ļ��д��ڣ�');
end
% % //========================================  
if isstr(fname)
   fname=cellstr(fname);
end
h=figure
c=colormap(jet(length(fname)+14));
sz = linspace(1,100,length(fname));
for i_file=1:length(fname)
Data = importdata(fullfile(location,char(fname(i_file)))); %ѡ���ļ���������
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
% figure;plot(real(cell2mat(omega)).');%˥������
% figure;plot(imag(cell2mat(omega)).');%Ƶ��--������ٶ�
axis equal
axis([-1.1 1.1 -1.1 1.1]);
