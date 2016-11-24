function devidBy3color(filename, filenum,startRow, endRow)
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
    figure(4*k-3)
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
    
    %% �O�F�\���i���x�j
    figure(4*k-2)
    L = mydata{k}.velocity <= 1.0;
    M = 1.0 < mydata{k}.velocity & mydata{k}.velocity <= 1.4; 
    N = 1.4 < mydata{k}.velocity;

    scatter3(mydata{k}.x(L),mydata{k}.y(L),mydata{k}.velocity(L),'b')
    hold on
    scatter3(mydata{k}.x(M),mydata{k}.y(M),mydata{k}.velocity(M),'g')
    scatter3(mydata{k}.x(N),mydata{k}.y(N),mydata{k}.velocity(N),'r')
    hold off

        %% �O�����O���t���i�Q�x�j
    figure(4*k-1)
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
    
    %% �O�F�\���i�Q�x�j
    figure(4*k)
    L = mydata{k}.vorticity <= 1.0;
    M = 1.0 < mydata{k}.velocity & mydata{k}.vorticity <= 1.4; 
    N = 1.4 < mydata{k}.vorticity;
    
    scatter3(mydata{k}.x(L),mydata{k}.y(L),mydata{k}.vorticity(L),'b')
    hold on
    scatter3(mydata{k}.x(M), mydata{k}.y(M), mydata{k}.vorticity(M),'g')
    scatter3(mydata{k}.x(N), mydata{k}.y(N), mydata{k}.vorticity(N),'r')
    hold off
    
end
