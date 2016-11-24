function T = dataCollectionGraph(filename, filenum,startRow, endRow)
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

for k = 1:numfiles
    %% 三次元グラフ化（速度）

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
    
        %% 三次元グラフ化（渦度）
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

%% 平均、分散、標準偏差、歪度、尖度（速度）

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
    v_mean(k) = mean(mydata{k}.velocity);    %平均
    v_var(k) = var(mydata{k}.velocity); %分散
    v_std(k) = std(mydata{k}.velocity); %標準偏差
    
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
vSkewness = v_skewness';
vKurtosis = v_kurtosis';

%% 平均、分散、標準偏差、歪度、尖度（渦度）

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
    v_mean(k) = mean(mydata{k}.vorticity);    %平均
    v_var(k) = var(mydata{k}.vorticity); %分散
    v_std(k) = std(mydata{k}.vorticity); %標準偏差
    
    % n = dataの長さ
    n = length(mydata{k}.vorticity);
    
    for q = 1:n
    s = s +(mydata{k}.vorticity(q) - v_mean(k))^3;
    v = v +(mydata{k}.vorticity(q) - v_mean(k))^4;
    end

    v_skewness(k) = s / (n * sqrt(v_var(k))^3); %歪度
    v_kurtosis(k) = v / (n * v_var(k)^2); %尖度
     
end

uMean = v_mean';
uVar = v_var';
uStd = v_std';
uSkewness = v_skewness';
uKurtosis = v_kurtosis';

%% Excelファイルに出力
%tabale作り
T = table(vMean, vVar, vStd, vSkewness, vKurtosis, uMean, uVar, uStd, uSkewness, uKurtosis, 'RowNames', fileName);
writetable(T, 'test.xls','WriteRowNames',true) %Row名付でファイル出力
