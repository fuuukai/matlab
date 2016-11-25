function T = dataCollectionGraph(filename, filenum,startRow, endRow)
%%
%IMPORTFILE1 �e�L�X�g �t�@�C�����琔�l�f�[�^���s��Ƃ��ăC���|�[�g���܂��B
%   F201POST1 = IMPORTFILE1(FILENAME) ����̑I���ɂ��Ă� �e�L�X�g �t�@�C�� FILENAME
%   ����f�[�^��ǂݎ��܂��B
%
%   F201POST1 = IMPORTFILE1(FILENAME, STARTROW, ENDROW) �e�L�X�g �t�@�C�� FILENAME
%   �� STARTROW �s���� ENDROW �s�܂ł̃f�[�^��ǂݎ��܂��B
%
% Example:
%   F201post1 = importfile1('F201post 1.csv', 2, 161);
%
%    TEXTSCAN ���Q�Ƃ��Ă��������B

% MATLAB �ɂ�鎩������ 2016/11/08 14:11:26

%% �ϐ������������܂��B
delimiter = ',';
if nargin<=2
    startRow = 2;
    endRow = inf;
end

%% �e�L�X�g�̊e�s�̏���������:
%   ��1: double (%f)
%	��2: double (%f)
%   ��3: double (%f)
%	��4: double (%f)
%   ��5: double (%f)
%	��6: double (%f)
%   ��7: double (%f)
%	��8: double (%f)
%   ��9: double (%f)
%	��10: double (%f)
%   ��11: double (%f)
%	��12: double (%f)
%   ��13: double (%f)
%	��14: �e�L�X�g (%s)
%   ��15: double (%f)
%	��16: double (%f)
% �ڍׂ� TEXTSCAN �̃h�L�������e�[�V�������Q�Ƃ��Ă��������B
formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%s%f%f%[^\n\r]';

%% �A�������t�@�C���̓ǂݍ���
%�@�t�@�C���̐����w�肷��
numfiles = filenum;
myfilename = cell(1,numfiles);
mydata = cell(1, numfiles);

for k = 1:numfiles
  myfilename{k} = sprintf('%s%d.csv',filename, k);

    %% �e�L�X�g �t�@�C�����J���܂��B
    fileID = fopen(myfilename{k},'r');

    % �f�[�^�̗������������ɏ]���ēǂݎ��܂��B
    % ���̌Ăяo���́A���̃R�[�h�̐����Ɏg�p���ꂽ�t�@�C���̍\���Ɋ�Â��Ă��܂��B�ʂ̃t�@�C���ŃG���[����������ꍇ�́A�C���|�[�g
    % �c�[������R�[�h�̍Đ��������݂Ă��������B
    dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
    for block=2:length(startRow)
        frewind(fileID);
        dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'HeaderLines', startRow(block)-1, 'ReturnOnError', false);
        for col=1:length(dataArray)
            dataArray{col} = [dataArray{col};dataArrayBlock{col}];
        end
    end
    
    %% �e�L�X�g �t�@�C������܂��B
    fclose(fileID);

    % �C���|�[�g�ł��Ȃ��f�[�^�̌㏈���B
    % �C���|�[�g���ɁA�C���|�[�g�ł��Ȃ��f�[�^�̋K�����K�p����Ȃ��������߁A�㏈���R�[�h���܂܂�Ă��܂���B�C���|�[�g�ł��Ȃ��f�[�^�ɓK�p�ł���R�[�h�𐶐�����ɂ́A�t�@�C�����̃C���|�[�g�ł��Ȃ��Z����I�����Ă���X�N���v�g���Đ������܂��B

    %% �o�͕ϐ��̍쐬
    mydata{k} = table(dataArray{1:end-1}, 'VariableNames', {'x','y','z','static_pressure','density','velocity','velocity_x','velocity_y','velocity_z','temparatuer','vorticity','total_pressuer','dinamic_pressure','share_stress','relative_pressure','viscosity'});
  
end

for k = 1:numfiles
    %% �O�����O���t���i���x�j

    figure(2*k-1)
    x = mydata{k}.x;                              % x axis
    y = mydata{k}.y;                                % y axis
    v = mydata{k}.velocity;                       % velocity

    scatter3(x,y,v,40,v,'filled')    % draw the scatter plot
    ax = gca;
    ax.XDir = 'reverse';
    view(-31,14)
    title('velocity on ' )
    xlabel('x axis [mm]')
    ylabel('y axis [mm]')
    zlabel('velocity[m/s]')

    cb = colorbar;                                     % create and label the colorbar
    cb.Label.String = 'veloicty [m/s]';
    
        %% �O�����O���t���i�Q�x�j
    figure(2*k)
    %x = fileName.x;                              % x axis
    %y = fileName.y;                                % y axis
    u = mydata{k}.vorticity;                       % vorticity

    scatter3(x,y,u,40,u,'filled')    % draw the scatter plot
    ax = gca;
    ax.XDir = 'reverse';
    view(-31,14)
    title('velocity on ' )
    xlabel('x axis [mm]')
    ylabel('y axis [mm]')
    zlabel('vorticity[1/s]')

    cb = colorbar;                                     % create and label the colorbar
    cb.Label.String = 'vorticity [1/s]';
    
end

%% ���ρA���U�A�W���΍��A�c�x�A��x�i���x�j

v_mean = zeros(1, numfiles);
v_var = zeros(1, numfiles);
v_std = zeros(1, numfiles);
fileName = {'a', 'b', 'c', 'd', 'e', 'f'};

s = 0;
v_skewness = zeros(1,numfiles); 
v = 0;
v_kurtosis = zeros(1, numfiles);
    
for k = 1:numfiles 
    fileName{k} = sprintf('%s%d.csv', filename, k);  
    v_mean(k) = mean(mydata{k}.velocity);    %����
    v_var(k) = var(mydata{k}.velocity); %���U
    v_std(k) = std(mydata{k}.velocity); %�W���΍�
    
    % n = data�̒���
    n = length(mydata{k}.velocity);
    
    for q = 1:n
    s = s +(mydata{k}.velocity(q) - v_mean(k))^3;
    v = v +(mydata{k}.velocity(q) - v_mean(k))^4;
    end

    v_skewness(k) = s / (n * sqrt(v_var(k))^3); %�c�x
    v_kurtosis(k) = v / (n * v_var(k)^2); %��x
     
end

vMean = v_mean';
vVar = v_var';
vStd = v_std';
vSkewness = v_skewness';
vKurtosis = v_kurtosis';

%% ���ρA���U�A�W���΍��A�c�x�A��x�i�Q�x�j

v_mean = zeros(1, numfiles);
v_var = zeros(1, numfiles);
v_std = zeros(1, numfiles);
fileName = {'a', 'b', 'c', 'd', 'e', 'f'};

s = 0;
v_skewness = zeros(1,numfiles); 
v = 0;
v_kurtosis = zeros(1, numfiles);
    
for k = 1:numfiles 
    fileName{k} = sprintf('%s%d.csv', filename, k);  
    v_mean(k) = mean(mydata{k}.vorticity);    %����
    v_var(k) = var(mydata{k}.vorticity); %���U
    v_std(k) = std(mydata{k}.vorticity); %�W���΍�
    
    % n = data�̒���
    n = length(mydata{k}.vorticity);
    
    for q = 1:n
    s = s +(mydata{k}.vorticity(q) - v_mean(k))^3;
    v = v +(mydata{k}.vorticity(q) - v_mean(k))^4;
    end

    v_skewness(k) = s / (n * sqrt(v_var(k))^3); %�c�x
    v_kurtosis(k) = v / (n * v_var(k)^2); %��x
     
end

uMean = v_mean';
uVar = v_var';
uStd = v_std';
uSkewness = v_skewness';
uKurtosis = v_kurtosis';

%% Excel�t�@�C���ɏo��
%tabale���
T = table(vMean, vVar, vStd, vSkewness, vKurtosis, uMean, uVar, uStd, uSkewness, uKurtosis, 'RowNames', fileName);
writetable(T, 'test.xls','WriteRowNames',true) %Row���t�Ńt�@�C���o��
