function T = dataCollectionNum(filename, filenum, startRow, endRow)
%%
%IMPORTFILE1 テキスト ファイルから数値データを行列としてインポートします。
%   F201POST1 = IMPORTFILE1(FILENAME) 既定の選択については テキスト ファイル FILENAME
%   からデータを読み取ります。
%
%   F201POST1 = IMPORTFILE1(FILENAME, STARTROW, ENDROW) テキスト ファイル FILENAME
%   の STARTROW 行から ENDROW 行までのデータを読み取ります。
%
% Example:
%   F201post1 = importfile1('F201post 1.csv', 2, 161);
%
%    TEXTSCAN も参照してください。

% MATLAB による自動生成 2016/11/08 14:11:26

%% 変数を初期化します。
delimiter = ',';
if nargin<=2
    startRow = 2;
    endRow = inf;
end

%% テキストの各行の書式文字列:
%   列1: double (%f)
%	列2: double (%f)
%   列3: double (%f)
%	列4: double (%f)
%   列5: double (%f)
%	列6: double (%f)
%   列7: double (%f)
%	列8: double (%f)
%   列9: double (%f)
%	列10: double (%f)
%   列11: double (%f)
%	列12: double (%f)
%   列13: double (%f)
%	列14: テキスト (%s)
%   列15: double (%f)
%	列16: double (%f)
% 詳細は TEXTSCAN のドキュメンテーションを参照してください。
formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%s%f%f%[^\n\r]';

%% 連続したファイルの読み込み
%　ファイルの数を指定する
numfiles = filenum;
myfilename = cell(1,numfiles);
mydata = cell(1, numfiles);

for k = 1:numfiles
  myfilename{k} = sprintf('%s%d.csv',filename, k);

    %% テキスト ファイルを開きます。
    fileID = fopen(myfilename{k},'r');

    % データの列を書式文字列に従って読み取ります。
    % この呼び出しは、このコードの生成に使用されたファイルの構造に基づいています。別のファイルでエラーが発生する場合は、インポート
    % ツールからコードの再生成を試みてください。
    dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
    for block=2:length(startRow)
        frewind(fileID);
        dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'HeaderLines', startRow(block)-1, 'ReturnOnError', false);
        for col=1:length(dataArray)
            dataArray{col} = [dataArray{col};dataArrayBlock{col}];
        end
    end
    
    %% テキスト ファイルを閉じます。
    fclose(fileID);

    % インポートできないデータの後処理。
    % インポート中に、インポートできないデータの規則が適用されなかったため、後処理コードが含まれていません。インポートできないデータに適用できるコードを生成するには、ファイル内のインポートできないセルを選択してからスクリプトを再生成します。

    %% 出力変数の作成
    mydata{k} = table(dataArray{1:end-1}, 'VariableNames', {'x','y','z','static_pressure','density','velocity','velocity_x','velocity_y','velocity_z','temparatuer','vorticity','total_pressuer','dinamic_pressure','share_stress','relative_pressure','viscosity'});
  
end


%% 平均、分散、標準偏差、変動係数、歪度、尖度（速度）

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
    v_mean(k) = mean(mydata{k}.velocity);    %平均
    v_var(k) = var(mydata{k}.velocity); %分散
    v_std(k) = std(mydata{k}.velocity); %標準偏差
    v_coefficientVariation(k) = v_std(k)/v_mean(k); %変動係数
    
    % n = dataの長さ
    n = length(mydata{k}.velocity);
    
    for q = 1:n
    s = s +(mydata{k}.velocity(q) - v_mean(k))^3;
    v = v +(mydata{k}.velocity(q) - v_mean(k))^4;
    end

    v_skewness(k) = s / (n * sqrt(v_var(k))^3); %歪度
    v_kurtosis(k) = v / (n * v_var(k)^2); %尖度
     
end

vMean = v_mean';
vVar = v_var';
vStd = v_std';
vCoefficientVariation = v_coefficientVariation';
vSkewness = v_skewness';
vKurtosis = v_kurtosis';

%% 平均、分散、標準偏差、変動係数、歪度、尖度（速度x）

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
    vx_mean(k) = mean(mydata{k}.velocity_x);    %平均
    vx_var(k) = var(mydata{k}.velocity_x); %分散
    vx_std(k) = std(mydata{k}.velocity_x); %標準偏差
    vx_coefficientVariation(k) = vx_std(k)/vx_mean(k); %変動係数
    
    % n = dataの長さ
    n = length(mydata{k}.velocity_x);
    
    for q = 1:n
    s = s +(mydata{k}.velocity_x(q) - vx_mean(k))^3;
    v = v +(mydata{k}.velocity_x(q) - vx_mean(k))^4;
    end

    vx_skewness(k) = s / (n * sqrt(vx_var(k))^3); %歪度
    vx_kurtosis(k) = v / (n * vx_var(k)^2); %尖度
     
end

vxMean = vx_mean';
vxVar = vx_var';
vxStd = vx_std';
vxCoefficientVariation = vx_coefficientVariation';
vxSkewness = vx_skewness';
vxKurtosis = vx_kurtosis';

%% 平均、分散、標準偏差、変動係数、歪度、尖度（速度y）

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
    vy_mean(k) = mean(mydata{k}.velocity_y);    %平均
    vy_var(k) = var(mydata{k}.velocity_y); %分散
    vy_std(k) = std(mydata{k}.velocity_y); %標準偏差
    vy_coefficientVariation(k) = vy_std(k)/vy_mean(k); %変動係数
    
    % n = dataの長さ
    n = length(mydata{k}.velocity_y);
    
    for q = 1:n
    s = s +(mydata{k}.velocity_y(q) - vy_mean(k))^3;
    v = v +(mydata{k}.velocity_y(q) - vy_mean(k))^4;
    end

    vy_skewness(k) = s / (n * sqrt(vy_var(k))^3); %歪度
    vy_kurtosis(k) = v / (n * vy_var(k)^2); %尖度
     
end

vyMean = vy_mean';
vyVar = vy_var';
vyStd = vy_std';
vyCoefficientVariation = vy_coefficientVariation';
vySkewness = vy_skewness';
vyKurtosis = vy_kurtosis';

%% 平均、分散、標準偏差、変動係数、歪度、尖度（速度z）

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
    vz_mean(k) = mean(mydata{k}.velocity_z);    %平均
    vz_var(k) = var(mydata{k}.velocity_z); %分散
    vz_std(k) = std(mydata{k}.velocity_z); %標準偏差
    vz_coefficientVariation(k) = vz_std(k)/vz_mean(k); %変動係数
    
    % n = dataの長さ
    n = length(mydata{k}.velocity_z);
    
    for q = 1:n
    s = s +(mydata{k}.velocity_z(q) - vz_mean(k))^3;
    v = v +(mydata{k}.velocity_z(q) - vz_mean(k))^4;
    end

    vz_skewness(k) = s / (n * sqrt(vz_var(k))^3); %歪度
    vz_kurtosis(k) = v / (n * vz_var(k)^2); %尖度
     
end

vzMean = vz_mean';
vzVar = vz_var';
vzStd = vz_std';
vzCoefficientVariation = vz_coefficientVariation';
vzSkewness = vz_skewness';
vzKurtosis = vz_kurtosis';

%% 平均、分散、標準偏差、変動係数、歪度、尖度（渦度）

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
    u_mean(k) = mean(mydata{k}.vorticity);    %平均
    u_var(k) = var(mydata{k}.vorticity); %分散
    u_std(k) = std(mydata{k}.vorticity); %標準偏差
    u_coefficientVariation(k) = u_std(k)/u_mean(k); %変動係数
    
    % n = dataの長さ
    n = length(mydata{k}.vorticity);
    
    for q = 1:n
    s = s +(mydata{k}.vorticity(q) - u_mean(k))^3;
    v = v +(mydata{k}.vorticity(q) - u_mean(k))^4;
    end

    u_skewness(k) = s / (n * sqrt(v_var(k))^3); %歪度
    u_kurtosis(k) = v / (n * v_var(k)^2); %尖度
     
end

uMean = u_mean';
uVar = u_var';
uStd = u_std';
uCoefficientVariation = u_coefficientVariation';
uSkewness = u_skewness';
uKurtosis = u_kurtosis';

%% 平均、分散、標準偏差、変動係数、歪度、尖度（全圧）

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
    tp_mean(k) = mean(mydata{k}.total_pressuer);    %平均
    tp_var(k) = var(mydata{k}.total_pressuer); %分散
    tp_std(k) = std(mydata{k}.total_pressuer); %標準偏差
    tp_coefficientVariation(k) = tp_std(k)/tp_mean(k); %変動係数
    
    % n = dataの長さ
    n = length(mydata{k}.total_pressuer);
    
    for q = 1:n
    s = s +(mydata{k}.total_pressuer(q) - tp_mean(k))^3;
    v = v +(mydata{k}.total_pressuer(q) - tp_mean(k))^4;
    end

    tp_skewness(k) = s / (n * sqrt(v_var(k))^3); %歪度
    tp_kurtosis(k) = v / (n * v_var(k)^2); %尖度
     
end

tpMean = tp_mean';
tpVar = tp_var';
tpStd = tp_std';
tpCoefficientVariation = tp_coefficientVariation';
tpSkewness = tp_skewness';
tpKurtosis = tp_kurtosis';

%% 平均、分散、標準偏差、変動係数、歪度、尖度（静圧）

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
    sp_mean(k) = mean(mydata{k}.static_pressure);    %平均
    sp_var(k) = var(mydata{k}.static_pressure); %分散
    sp_std(k) = std(mydata{k}.static_pressure); %標準偏差
    sp_coefficientVariation(k) = sp_std(k)/sp_mean(k); %変動係数
    
    % n = dataの長さ
    n = length(mydata{k}.static_pressure);
    
    for q = 1:n
    s = s +(mydata{k}.static_pressure(q) - sp_mean(k))^3;
    v = v +(mydata{k}.static_pressure(q) - sp_mean(k))^4;
    end

    sp_skewness(k) = s / (n * sqrt(v_var(k))^3); %歪度
    sp_kurtosis(k) = v / (n * v_var(k)^2); %尖度
     
end

spMean = sp_mean';
spVar = sp_var';
spStd = sp_std';
spCoefficientVariation = sp_coefficientVariation';
spSkewness = sp_skewness';
spKurtosis = sp_kurtosis';

%% 平均、分散、標準偏差、変動係数、歪度、尖度（動圧）

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
    dp_mean(k) = mean(mydata{k}.dinamic_pressure);    %平均
    dp_var(k) = var(mydata{k}.dinamic_pressure); %分散
    dp_std(k) = std(mydata{k}.dinamic_pressure); %標準偏差
    dp_coefficientVariation(k) = dp_std(k)/dp_mean(k); %変動係数
    
    % n = dataの長さ
    n = length(mydata{k}.dinamic_pressure);
    
    for q = 1:n
    s = s +(mydata{k}.dinamic_pressure(q) - dp_mean(k))^3;
    v = v +(mydata{k}.dinamic_pressure(q) - dp_mean(k))^4;
    end

    dp_skewness(k) = s / (n * sqrt(v_var(k))^3); %歪度
    dp_kurtosis(k) = v / (n * v_var(k)^2); %尖度
     
end

dpMean = dp_mean';
dpVar = dp_var';
dpStd = dp_std';
dpCoefficientVariation = dp_coefficientVariation';
dpSkewness = dp_skewness';
dpKurtosis = dp_kurtosis';

%% Excelファイルに出力
%tabale作り
T = table(vMean, vVar, vStd, vCoefficientVariation, vSkewness, vKurtosis, vxMean, vxVar, vxStd, vxCoefficientVariation, vxSkewness, vxKurtosis, vyMean, vyVar, vyStd, vyCoefficientVariation, vySkewness, vyKurtosis, vzMean, vzVar, vzStd, vzCoefficientVariation, vzSkewness, vzKurtosis, tpMean, tpVar, tpStd, tpCoefficientVariation, tpSkewness, tpKurtosis, spMean, spVar, spStd, spCoefficientVariation, spSkewness, spKurtosis, dpMean, dpVar, dpStd, dpCoefficientVariation, dpSkewness, dpKurtosis, uMean, uVar, uStd, uCoefficientVariation, uSkewness, uKurtosis, 'RowNames', fileName);
writetable(T, 'test.xls','WriteRowNames',true) %Row名付でファイル出力
