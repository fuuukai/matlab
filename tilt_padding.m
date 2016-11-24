%% 読み込み
function tableMake = tilt_padding(filename, paramaterColumnNum, Dx, Dy)

%% テーブルの作成

csvFileName = sprintf('%s.csv', filename);
T = readtable(csvFileName);

%% x,yの重複なしで要素抽出

Tx = round(T.x,5);
Ty = round(T.y,5);

xComponent = unique(double(Tx));
yComponent = unique(double(Ty));
disp('xComponent = ')
disp(xComponent)
disp('yComponent = ')
disp(yComponent)
%%

xCompoArray = double(xComponent');
yCompoArray = double(yComponent');
disp('xCompoArray = ')
disp(xCompoArray)
disp('yCompoArray = ')
disp(yCompoArray)
%% 配列の作成と並び替え

xLength = length(xCompoArray);
%disp(xLength)
yLength = length(yCompoArray);
%disp(yLength)
dataNum = length(Tx);

%% 探索と新たな配列作成ベクトルデータ格納

image = cell(yLength, xLength);
disp(image)

for c = 1:dataNum
    
    for r = 1:yLength
        if yCompoArray(1,r) == round(T.y(c), 5)
            tmpY = r;
            break
        end
    end
    
    for q = 1:xLength
        if xCompoArray(1,q) == round(T.x(c), 5)
            tmpX = q;
            break
        end
    end
    
    image{tmpY, tmpX} = T(c,:);

end
    disp('image = ')
    disp(image)
    

%% １つのパラメータ抽出

%nameParamater = char('static_pressure', 'velocity', 'velocity_x', 'velocity_y', 'velocity_z', 'vorticity', 'total_pressure', 'dinamic_pressure');
%numParamater = [4, 6, 7, 8, 9, 11, 12, 13];

imV  = zeros(yLength+2, xLength+2);

for s = 2:yLength+1
    for t = 2:xLength+1
        if ~isempty(image{s-1,t-1})
            disp('s=')
            disp(s)
            disp('t=')
            disp(t)
            tmp = table2array((image{s-1,t-1}) );
            imV(s,t) = tmp(1, paramaterColumnNum);

        elseif isempty(image{s-1,t-1})
            imV(s,t) = 0;
        end
    end
end
    time = cputime;
    disp(time)
disp(imV)


%% 傾斜量

tiltMap = zeros(yLength, xLength);

for j = 2:yLength+1
    for i = 2:xLength+1
        
        Sx = ( imV(j-1, i-1) + imV(j,i-1) + imV(j+1,i-1) - imV(j-1,i+1) - imV(j,i+1) - imV(j+1,i+1) ) / 6 * Dx;
        Sy = ( imV(j-1,i-1) + imV(j-1,i) + imV(j-1,i+1) - imV(j+1, i-1) - imV(j+1, i) - imV(j+1,i+1) ) / 6 * Dy;
        
        tiltMap(j-1,i-1) = sqrt(Sx^2 + Sy^2);

    end
end
        disp(tiltMap)
        imshow(tiltMap)

%% 斜面曲率

%斜面曲率
slopeCurvature = zeros(yLength, xLength);
%曲率
curvature = zeros(yLength, xLength);

for j = 2:yLength+1
    for i = 2:xLength+1
        % A = [(Z1 + Z3 + Z7 + Z9) / 4  - (Z2 + Z4 + Z6 + Z8) / 2 + Z5]/ L^4 
        % B = [(Z1 + Z3 - Z7 - Z9) /4 - (Z2 - Z8) /2] / L^3
        % C = [(-Z1 + Z3 - Z7 + Z9) /4 + (Z4 - Z6)] /2] / L^3
        % D = [(Z4 + Z6) /2 - Z5] / L^2
        % E = [(Z2 + Z8) /2 - Z5] / L2
        % F = (-Z1 + Z3 + Z7 - Z9) / 4L^2 
        % G = (-Z4 + Z6) / 2L
        % H = (Z2 - Z8) / 2L
        % II = Z5
       
        %yとxの座標反対になっている気がする
        %a = ( ( t(i-1,j-1) + t(i+1, j-1) + t(i-1, j+1) + t(i+1, j+1) ) / 4 - ( t() + t() + t() + t()) / 2 + t(i, j) ) / Dx^4;
        %b = ( ( t(i-1,j-1) + t(i+1, j-1) - t(i-1, j+1) - t(i+1, j+1) ) / 4 - ( t(i, j-1) + t(i, j+1) / 2 ) ) / Dx^3;
        %c = ( ( -t(i-1,j-1) + t(i+1, j-1) - t(i-1, j+1) + t(i+1, j+1) ) / 4 - ( t(i-1, j) - t(i+1, j) / 2 ) ) / Dx^3;
        %d = ( ( t(i-1, j) + t(i+1, j) )/ 2 - t(i,j) )  / Dx^2;
        %e = ( ( t(i, j-1) + t(i, j+1) )/ 2 - t(i,j) )  / Dx^2;
        %f = ( ( -t(i-1, j-1) + t(i+1, j-1) ) + t(i-1, j+1) - t(i+1, j+1) )  / 4*Dx^2;
        %g = ( -t(i-1, j) + t(i+1, j) ) * 2*Dx;
        %h = ( t(i, j-1) - t(i, j+1) ) * 2*Dx;
        %ii = t(i, j);
        
        a = ( ( imV(j-1, i-1) + imV(j-1, i+1) + imV(j+1, i-1) + imV(j+1, i+1) ) / 4 - ( imV(j-1, i) + imV(j, i-1) + imV(j, i+1) + imV(j+1, i)) / 2 + imV(j, i) ) / Dx^4;
        b = ( ( imV(j-1,i-1) + imV(j-1, i+1) - imV(j+1, i-1) - imV(j+1, i+1) ) / 4 - ( imV(j-1, i) + imV(j+1, i) / 2 ) ) / Dx^3;
        c = ( ( -imV(j-1,i-1) + imV(j-1, i+1) - imV(j+1, i-1) + imV(j+1, i+1) ) / 4 - ( imV(j, i-1) - imV(j, i+1) / 2 ) ) / Dx^3;       
        d = ( ( imV(j,i-1) + imV(j, i+1) )/ 2 - imV(j,i) )  / Dx^2;
        e = ( ( imV(j-1,i) + imV(j+1, i) )/ 2 - imV(j,i) )  / Dx^2;
        
        f = ( ( -imV(j-1,i-1) + imV(j-1,i+1) ) + imV(j+1,i-1) - imV(j+1, i+1) )  / 4*Dx^2;
        g = ( -imV(j,i-1) + imV(j,i+1) ) * 2*Dx;
        h = ( imV(j-1,i) - imV(j+1,i) ) * 2*Dx;
        ii = imV(j,i);
        
        slopeCurvature = 0;
        curvature(j-1,i-1) = -2 * (d+e) * 100;
        
    end
end

%% 斜面方位図

slopeOrientation = zeros(yLength, xLength);

for j = 2:yLength+1
    for i = 2:xLength+1

        Hx = ( imV(j-1,i-1) + imV(j,i-1) + imV(j+1,i-1) - ( imV(j-1, i+1) + imV(j,i+1) + imV(j+1,i+1) ) ) / (3*Dx);
        Hy = ( imV(j+1,i-1) + imV(j+1,i) + imV(j+1,i+1) - ( imV(j-1,i-1) + imV(j-1,i) + imV(j-1,i+1) ) ) / (3*Dx);
        
        slopeOrientation(j-1,i-1) =  atan( (Hx/Hy)*180/pi );
        
    end
end
disp('slopeOrientation=')
disp(slopeOrientation)
%% テーブル作成

tableMake = table(tiltMap, curvature, slopeOrientation);
nameTiltFile = sprintf('%s_tilt.xls', filename);
writetable(tableMake, nameTiltFile)

%% グラフ表示

figure(1)
imshow(tableMake.tiltMap)
figure(2)
histogram(tableMake.tiltMap)
figure(3)
histogram(tableMake.curvature)
figure(4)
histogram(tableMake.slopeOrientation)

end
