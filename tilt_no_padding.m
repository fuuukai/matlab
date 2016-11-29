%% �ǂݍ���
function tableMake = tilt_no_padding(filename, paramaterName, paramaterColumnNum, Dx, Dy)

%% �e�[�u���̍쐬

csvFileName = sprintf('%s.csv', filename);
T = readtable(csvFileName);

%% x,y�̏d���Ȃ��ŗv�f���o

Tx = round(T.x,5);
Ty = round(T.y,5);

arrayLength = length(Tx);

xComponent = unique(double(Tx));
yComponent = unique(double(Ty));
disp('xComponent = ')
disp(xComponent)
disp('yComponent = ')
disp(yComponent)

xCompoArray = double(xComponent');
yCompoArray = double(yComponent');
disp('xCompoArray = ')
disp(xCompoArray)
disp('yCompoArray = ')
disp(yCompoArray)
%% �z��̍쐬�ƕ��ёւ�

xLength = length(xCompoArray);
%disp(xLength)
yLength = length(yCompoArray);
%disp(yLength)
dataNum = length(Tx);

%% �T���ƐV���Ȕz��쐬�x�N�g���f�[�^�i�[

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
    

%% �P�̃p�����[�^���o

%nameParamater = char('static_pressure', 'velocity', 'velocity_x', 'velocity_y', 'velocity_z', 'vorticity', 'total_pressure', 'dinamic_pressure');
%numParamater = [4, 6, 7, 8, 9, 11, 12, 13];

imV  = zeros(yLength, xLength);

for s = 1:yLength
    for t = 1:xLength
        if ~isempty(image{s,t})
            disp('s=')
            disp(s)
            disp('t=')
            disp(t)
            tmp = table2array((image{s,t}) );
            imV(s,t) = tmp(1, paramaterColumnNum);

        elseif isempty(image{s,t})
            imV(s,t) = 0;
        end
    end
end
    time = cputime;
    disp(time)
disp(imV)
figure(1)


%% �X�Η�

tiltMap = zeros(yLength, xLength);
tiltMapArray = NaN(arrayLength, 1);
count =1;

for j = 2:yLength-1
    for i = 2:xLength-1
        
        Sx = ( imV(j-1, i-1) + imV(j,i-1) + imV(j+1,i-1) - imV(j-1,i+1) - imV(j,i+1) - imV(j+1,i+1) ) / 6 * Dx;
        Sy = ( imV(j-1,i-1) + imV(j-1,i) + imV(j-1,i+1) - imV(j+1, i-1) - imV(j+1, i) - imV(j+1,i+1) ) / 6 * Dy;
        
        tiltMap(j,i) = sqrt(Sx^2 + Sy^2);
        
        if imV(j,i) ~= 0
            tiltMapArray(count, 1) = tiltMap(j,i);
            count = count+1;
        end
            
    end
end

b = ~isnan(tiltMapArray);
tiltMapArray = tiltMapArray(b);

disp(tiltMap)
imshow(tiltMap)
disp(tiltMapArray)

%% �Ζʋȗ�

%�Ζʋȗ�
slopeCurvature = zeros(yLength, xLength);
%�ȗ�
curvature = zeros(yLength, xLength);

curvatureArray = NaN(arrayLength, 1);
count =1;

for j = 2:yLength-1
    for i = 2:xLength-1
        % A = [(Z1 + Z3 + Z7 + Z9) / 4  - (Z2 + Z4 + Z6 + Z8) / 2 + Z5]/ L^4 
        % B = [(Z1 + Z3 - Z7 - Z9) /4 - (Z2 - Z8) /2] / L^3
        % C = [(-Z1 + Z3 - Z7 + Z9) /4 + (Z4 - Z6)] /2] / L^3
        % D = [(Z4 + Z6) /2 - Z5] / L^2
        % E = [(Z2 + Z8) /2 - Z5] / L2
        % F = (-Z1 + Z3 + Z7 - Z9) / 4L^2 
        % G = (-Z4 + Z6) / 2L
        % H = (Z2 - Z8) / 2L
        % II = Z5
       
        %y��x�̍��W���΂ɂȂ��Ă���C������
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
        curvature(j,i) = -2 * (d+e) * 100;
        
        if imV(j,i) ~= 0
            curvatureArray(count, 1) = curvature(j,i);
            count = count+1;
        end
     
    end
end

b = ~isnan(curvatureArray);
curvatureArray = curvatureArray(b);
disp(curvatureArray)

%% �Ζʕ��ʐ}

slopeOrientation = zeros(yLength, xLength);

slopeOrientationArray = NaN(arrayLength, 1);
count =1;

for j = 2:yLength-1
    for i = 2:xLength-1

        Hx = ( imV(j-1,i-1) + imV(j,i-1) + imV(j+1,i-1) - ( imV(j-1, i+1) + imV(j,i+1) + imV(j+1,i+1) ) ) / (3*Dx);
        Hy = ( imV(j+1,i-1) + imV(j+1,i) + imV(j+1,i+1) - ( imV(j-1,i-1) + imV(j-1,i) + imV(j-1,i+1) ) ) / (3*Dx);
        
        slopeOrientation(j,i) =  atan( (Hx/Hy)*180/pi );
        
        if imV(j,i) ~= 0
            slopeOrientationArray(count, 1) = slopeOrientation(j,i);
            count = count+1;
        end
        
    end
end

b = ~isnan(slopeOrientationArray);
slopeOrientationArray = slopeOrientationArray(b);
disp(slopeOrientationArray)

%% �e�[�u���쐬

tableMake = table(tiltMap, curvature, slopeOrientation);
nameTiltFile = sprintf('%s_tilt_no_padding.xlsx', filename);
writetable(tableMake, nameTiltFile)

tableMakeNoNanData = table(tiltMapArray, curvatureArray, slopeOrientationArray);
disp(tableMakeNoNanData)

%% �O���t�\���Ɛ}�̕ۑ�

paramatername = paramaterName;
nbins = 100;

figure(1)
imshow(tableMake.tiltMap)
title1 = sprintf('%s %s image', filename,paramatername);
figureName1 = sprintf('%s_%s_image_no_padding.png', filename,paramatername);
title(title1)
saveas(gcf, figureName1);

% figure(2)
% histogram(tableMake.tiltMap, nbins)
% title2 = sprintf('%s %s tiltMap', filename,paramatername);
% figureName2 = sprintf('%s_%s_tiltMap_no_padding.png', filename,paramatername);
% title(title2)
%saveas(gcf, figureName2);
% 
% figure(3)
% histogram(tableMake.curvature, nbins)
% title3 = sprintf('%s %s curvature', filename,paramatername);
% figureName3 = sprintf('%s_%s_curvature_no_padding.png', filename,paramatername);
% title(title3)
% saveas(gcf, figureName3);
% 
% figure(4)
% histogram(tableMake.slopeOrientation, nbins)
% title4 = sprintf('%s %s slopeOrientation', filename,paramatername);
% figureName4 = sprintf('%s_%s_slopeOrientation_no_padding.png', filename,paramatername);
% title(title4)
% saveas(gcf, figureName4);

figure(2)
histogram(tableMakeNoNanData.tiltMapArray, nbins)
title2 = sprintf('%s %s tiltMap', filename,paramatername);
figureName2 = sprintf('%s_%s_tiltMap_no_padding.png', filename,paramatername);
title(title2)
saveas(gcf, figureName2);

figure(3)
histogram(tableMakeNoNanData.curvatureArray, nbins)
title3 = sprintf('%s %s curvature', filename,paramatername);
figureName3 = sprintf('%s_%s_curvature_no_padding.png', filename,paramatername);
title(title3)
saveas(gcf, figureName3);


figure(4)
histogram(tableMakeNoNanData.slopeOrientationArray, nbins)
title5 = sprintf('%s %s slopeOrientationArray', filename,paramatername);
figureName5 = sprintf('%s_%s_slopeOrientation_no_padding.png', filename,paramatername);
title(title5)
saveas(gcf, figureName5);

%% ���ρA���U�A�W���΍��A�ϓ��W���A�c�x�A��x�i�X�Γx�j

til_mean = mean(tiltMapArray);    %����
til_var = var(tiltMapArray); %���U
til_std = std(tiltMapArray); %�W���΍�
til_coefficientVariation = til_std/til_mean; %�ϓ��W��

s=0;
v=0;

% n = data�̒���
n = length(tiltMapArray);

for q = 1:n
s = s +(tiltMapArray(q) - til_mean)^3;
v = v +(tiltMapArray(q) - til_mean)^4;
end

til_skewness = s / (n * sqrt(til_var)^3); %�c�x
til_kurtosis = v / (n * til_var^2); %��x


tilMean = til_mean';
tilVar = til_var';
tilStd = til_std';
tilCoefficientVariation = til_coefficientVariation';
tilSkewness = til_skewness';
tilKurtosis = til_kurtosis';

%% ���ρA���U�A�W���΍��A�ϓ��W���A�c�x�A��x�i�Ζʋȗ��j

cur_mean = mean(curvatureArray);    %����
cur_var = var(curvatureArray); %���U
cur_std = std(curvatureArray); %�W���΍�
cur_coefficientVariation = cur_std/cur_mean; %�ϓ��W��

s=0;
v=0;

% n = data�̒���
n = length(curvatureArray);

for q = 1:n
s = s +(curvatureArray(q) - cur_mean)^3;
v = v +(curvatureArray(q) - cur_mean)^4;
end

cur_skewness = s / (n * sqrt(cur_var)^3); %�c�x
cur_kurtosis = v / (n * cur_var^2); %��x


curMean = cur_mean';
curVar = cur_var';
curStd = cur_std';
curCoefficientVariation = cur_coefficientVariation';
curSkewness = cur_skewness';
curKurtosis = cur_kurtosis';

%% ���ρA���U�A�W���΍��A�ϓ��W���A�c�x�A��x�i�Ζʕ��ʐ}�j

ori_mean = mean(slopeOrientationArray);    %����
ori_var = var(slopeOrientationArray); %���U
ori_std = std(slopeOrientationArray); %�W���΍�
ori_coefficientVariation = ori_std/ori_mean; %�ϓ��W��

s=0;
v=0;

% n = data�̒���
n = length(slopeOrientationArray);

for q = 1:n
s = s +(slopeOrientationArray(q) - ori_mean)^3;
v = v +(slopeOrientationArray(q) - ori_mean)^4;
ends

ori_skewness = s / (n * sqrt(ori_var)^3); %�c�x
ori_kurtosis = v / (n * ori_var^2); %��x


oriMean = ori_mean';
oriVar = ori_var';
oriStd = ori_std';
oriCoefficientVariation = ori_coefficientVariation';
oriSkewness = ori_skewness';
oriKurtosis = ori_kurtosis';

%% Table�쐬
T = table(tilMean, tilVar, tilStd, tilCoefficientVariation, tilSkewness, tilKurtosis, curMean, curVar, curStd, curCoefficientVariation, curSkewness, curKurtosis, oriMean ,oriVar, oriStd, oriCoefficientVariation, oriSkewness, oriKurtosis);
disp(T)
end
