function T = dataCollectionNum(filename, filenum, startRow, endRow)
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


%% ���ρA���U�A�W���΍��A�ϓ��W���A�c�x�A��x�i���x�j

v_mean = zeros(1, numfiles);
v_var = zeros(1, numfiles);
v_std = zeros(1, numfiles);
v_coefficientVariation = zeros(1, numfiles);
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
    v_coefficientVariation(k) = v_std(k)/v_mean(k); %�ϓ��W��
    
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
vCoefficientVariation = v_coefficientVariation';
vSkewness = v_skewness';
vKurtosis = v_kurtosis';

%% ���ρA���U�A�W���΍��A�ϓ��W���A�c�x�A��x�i���xx�j

vx_mean = zeros(1, numfiles);
vx_var = zeros(1, numfiles);
vx_std = zeros(1, numfiles);
vx_coefficientVariation = zeros(1, numfiles);
fileName = {'a', 'b', 'c', 'd', 'e', 'f'};

s = 0;
vx_skewness = zeros(1,numfiles); 
v = 0;
vx_kurtosis = zeros(1, numfiles);
    
for k = 1:numfiles 
    fileName{k} = sprintf('%s%d.csv', filename, k);  
    vx_mean(k) = mean(mydata{k}.velocity_x);    %����
    vx_var(k) = var(mydata{k}.velocity_x); %���U
    vx_std(k) = std(mydata{k}.velocity_x); %�W���΍�
    vx_coefficientVariation(k) = vx_std(k)/vx_mean(k); %�ϓ��W��
    
    % n = data�̒���
    n = length(mydata{k}.velocity_x);
    
    for q = 1:n
    s = s +(mydata{k}.velocity_x(q) - vx_mean(k))^3;
    v = v +(mydata{k}.velocity_x(q) - vx_mean(k))^4;
    end

    vx_skewness(k) = s / (n * sqrt(vx_var(k))^3); %�c�x
    vx_kurtosis(k) = v / (n * vx_var(k)^2); %��x
     
end

vxMean = vx_mean';
vxVar = vx_var';
vxStd = vx_std';
vxCoefficientVariation = vx_coefficientVariation';
vxSkewness = vx_skewness';
vxKurtosis = vx_kurtosis';

%% ���ρA���U�A�W���΍��A�ϓ��W���A�c�x�A��x�i���xy�j

vy_mean = zeros(1, numfiles);
vy_var = zeros(1, numfiles);
vy_std = zeros(1, numfiles);
vy_coefficientVariation = zeros(1, numfiles);
fileName = {'a', 'b', 'c', 'd', 'e', 'f'};

s = 0;
vy_skewness = zeros(1,numfiles); 
v = 0;
vy_kurtosis = zeros(1, numfiles);
    
for k = 1:numfiles 
    fileName{k} = sprintf('%s%d.csv', filename, k);  
    vy_mean(k) = mean(mydata{k}.velocity_y);    %����
    vy_var(k) = var(mydata{k}.velocity_y); %���U
    vy_std(k) = std(mydata{k}.velocity_y); %�W���΍�
    vy_coefficientVariation(k) = vy_std(k)/vy_mean(k); %�ϓ��W��
    
    % n = data�̒���
    n = length(mydata{k}.velocity_y);
    
    for q = 1:n
    s = s +(mydata{k}.velocity_y(q) - vy_mean(k))^3;
    v = v +(mydata{k}.velocity_y(q) - vy_mean(k))^4;
    end

    vy_skewness(k) = s / (n * sqrt(vy_var(k))^3); %�c�x
    vy_kurtosis(k) = v / (n * vy_var(k)^2); %��x
     
end

vyMean = vy_mean';
vyVar = vy_var';
vyStd = vy_std';
vyCoefficientVariation = vy_coefficientVariation';
vySkewness = vy_skewness';
vyKurtosis = vy_kurtosis';

%% ���ρA���U�A�W���΍��A�ϓ��W���A�c�x�A��x�i���xz�j

vz_mean = zeros(1, numfiles);
vz_var = zeros(1, numfiles);
vz_std = zeros(1, numfiles);
vz_coefficientVariation = zeros(1, numfiles);
fileName = {'a', 'b', 'c', 'd', 'e', 'f'};

s = 0;
vz_skewness = zeros(1,numfiles); 
v = 0;
vz_kurtosis = zeros(1, numfiles);
    
for k = 1:numfiles 
    fileName{k} = sprintf('%s%d.csv', filename, k);  
    vz_mean(k) = mean(mydata{k}.velocity_z);    %����
    vz_var(k) = var(mydata{k}.velocity_z); %���U
    vz_std(k) = std(mydata{k}.velocity_z); %�W���΍�
    vz_coefficientVariation(k) = vz_std(k)/vz_mean(k); %�ϓ��W��
    
    % n = data�̒���
    n = length(mydata{k}.velocity_z);
    
    for q = 1:n
    s = s +(mydata{k}.velocity_z(q) - vz_mean(k))^3;
    v = v +(mydata{k}.velocity_z(q) - vz_mean(k))^4;
    end

    vz_skewness(k) = s / (n * sqrt(vz_var(k))^3); %�c�x
    vz_kurtosis(k) = v / (n * vz_var(k)^2); %��x
     
end

vzMean = vz_mean';
vzVar = vz_var';
vzStd = vz_std';
vzCoefficientVariation = vz_coefficientVariation';
vzSkewness = vz_skewness';
vzKurtosis = vz_kurtosis';

%% ���ρA���U�A�W���΍��A�ϓ��W���A�c�x�A��x�i�Q�x�j

u_mean = zeros(1, numfiles);
u_var = zeros(1, numfiles);
u_std = zeros(1, numfiles);
u_coefficientVariation = zeros(1, numfiles);
fileName = {'a', 'b', 'c', 'd', 'e', 'f'};

s = 0;
u_skewness = zeros(1,numfiles); 
v = 0;
u_kurtosis = zeros(1, numfiles);
    
for k = 1:numfiles 
    fileName{k} = sprintf('%s%d.csv', filename, k);  
    u_mean(k) = mean(mydata{k}.vorticity);    %����
    u_var(k) = var(mydata{k}.vorticity); %���U
    u_std(k) = std(mydata{k}.vorticity); %�W���΍�
    u_coefficientVariation(k) = u_std(k)/u_mean(k); %�ϓ��W��
    
    % n = data�̒���
    n = length(mydata{k}.vorticity);
    
    for q = 1:n
    s = s +(mydata{k}.vorticity(q) - u_mean(k))^3;
    v = v +(mydata{k}.vorticity(q) - u_mean(k))^4;
    end

    u_skewness(k) = s / (n * sqrt(v_var(k))^3); %�c�x
    u_kurtosis(k) = v / (n * v_var(k)^2); %��x
     
end

uMean = u_mean';
uVar = u_var';
uStd = u_std';
uCoefficientVariation = u_coefficientVariation';
uSkewness = u_skewness';
uKurtosis = u_kurtosis';

%% ���ρA���U�A�W���΍��A�ϓ��W���A�c�x�A��x�i�S���j

tp_mean = zeros(1, numfiles);
tp_var = zeros(1, numfiles);
tp_std = zeros(1, numfiles);
tp_coefficientVariation = zeros(1, numfiles);
fileName = {'a', 'b', 'c', 'd', 'e', 'f'};

s = 0;
tp_skewness = zeros(1,numfiles); 
v = 0;
tp_kurtosis = zeros(1, numfiles);
    
for k = 1:numfiles 
    fileName{k} = sprintf('%s%d.csv', filename, k);  
    tp_mean(k) = mean(mydata{k}.total_pressuer);    %����
    tp_var(k) = var(mydata{k}.total_pressuer); %���U
    tp_std(k) = std(mydata{k}.total_pressuer); %�W���΍�
    tp_coefficientVariation(k) = tp_std(k)/tp_mean(k); %�ϓ��W��
    
    % n = data�̒���
    n = length(mydata{k}.total_pressuer);
    
    for q = 1:n
    s = s +(mydata{k}.total_pressuer(q) - tp_mean(k))^3;
    v = v +(mydata{k}.total_pressuer(q) - tp_mean(k))^4;
    end

    tp_skewness(k) = s / (n * sqrt(v_var(k))^3); %�c�x
    tp_kurtosis(k) = v / (n * v_var(k)^2); %��x
     
end

tpMean = tp_mean';
tpVar = tp_var';
tpStd = tp_std';
tpCoefficientVariation = tp_coefficientVariation';
tpSkewness = tp_skewness';
tpKurtosis = tp_kurtosis';

%% ���ρA���U�A�W���΍��A�ϓ��W���A�c�x�A��x�i�È��j

sp_mean = zeros(1, numfiles);
sp_var = zeros(1, numfiles);
sp_std = zeros(1, numfiles);
sp_coefficientVariation = zeros(1, numfiles);
fileName = {'a', 'b', 'c', 'd', 'e', 'f'};

s = 0;
sp_skewness = zeros(1,numfiles); 
v = 0;
sp_kurtosis = zeros(1, numfiles);
    
for k = 1:numfiles 
    fileName{k} = sprintf('%s%d.csv', filename, k);  
    sp_mean(k) = mean(mydata{k}.static_pressure);    %����
    sp_var(k) = var(mydata{k}.static_pressure); %���U
    sp_std(k) = std(mydata{k}.static_pressure); %�W���΍�
    sp_coefficientVariation(k) = sp_std(k)/sp_mean(k); %�ϓ��W��
    
    % n = data�̒���
    n = length(mydata{k}.static_pressure);
    
    for q = 1:n
    s = s +(mydata{k}.static_pressure(q) - sp_mean(k))^3;
    v = v +(mydata{k}.static_pressure(q) - sp_mean(k))^4;
    end

    sp_skewness(k) = s / (n * sqrt(v_var(k))^3); %�c�x
    sp_kurtosis(k) = v / (n * v_var(k)^2); %��x
     
end

spMean = sp_mean';
spVar = sp_var';
spStd = sp_std';
spCoefficientVariation = sp_coefficientVariation';
spSkewness = sp_skewness';
spKurtosis = sp_kurtosis';

%% ���ρA���U�A�W���΍��A�ϓ��W���A�c�x�A��x�i�����j

dp_mean = zeros(1, numfiles);
dp_var = zeros(1, numfiles);
dp_std = zeros(1, numfiles);
dp_coefficientVariation = zeros(1, numfiles);
fileName = {'a', 'b', 'c', 'd', 'e', 'f'};

s = 0;
dp_skewness = zeros(1,numfiles); 
v = 0;
dp_kurtosis = zeros(1, numfiles);
    
for k = 1:numfiles 
    fileName{k} = sprintf('%s%d.csv', filename, k);  
    dp_mean(k) = mean(mydata{k}.dinamic_pressure);    %����
    dp_var(k) = var(mydata{k}.dinamic_pressure); %���U
    dp_std(k) = std(mydata{k}.dinamic_pressure); %�W���΍�
    dp_coefficientVariation(k) = dp_std(k)/dp_mean(k); %�ϓ��W��
    
    % n = data�̒���
    n = length(mydata{k}.dinamic_pressure);
    
    for q = 1:n
    s = s +(mydata{k}.dinamic_pressure(q) - dp_mean(k))^3;
    v = v +(mydata{k}.dinamic_pressure(q) - dp_mean(k))^4;
    end

    dp_skewness(k) = s / (n * sqrt(v_var(k))^3); %�c�x
    dp_kurtosis(k) = v / (n * v_var(k)^2); %��x
     
end

dpMean = dp_mean';
dpVar = dp_var';
dpStd = dp_std';
dpCoefficientVariation = dp_coefficientVariation';
dpSkewness = dp_skewness';
dpKurtosis = dp_kurtosis';

%% Excel�t�@�C���ɏo��
%tabale���
T = table(vMean, vVar, vStd, vCoefficientVariation, vSkewness, vKurtosis, vxMean, vxVar, vxStd, vxCoefficientVariation, vxSkewness, vxKurtosis, vyMean, vyVar, vyStd, vyCoefficientVariation, vySkewness, vyKurtosis, vzMean, vzVar, vzStd, vzCoefficientVariation, vzSkewness, vzKurtosis, tpMean, tpVar, tpStd, tpCoefficientVariation, tpSkewness, tpKurtosis, spMean, spVar, spStd, spCoefficientVariation, spSkewness, spKurtosis, dpMean, dpVar, dpStd, dpCoefficientVariation, dpSkewness, dpKurtosis, uMean, uVar, uStd, uCoefficientVariation, uSkewness, uKurtosis, 'RowNames', fileName);
writetable(T, 'test.xls','WriteRowNames',true) %Row���t�Ńt�@�C���o��
