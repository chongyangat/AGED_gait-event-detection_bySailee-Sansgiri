%% Main_Algorithms - Modified for Single C3D File Analysis (No Force Plate)
% Modified to analyze single C3D file: BWA7.c3d
% Adapted for available marker set without SACR marker

% Codes connected to Visscher & Sansgiri et al. "Towards validation and 
% standardization of automatic gait event identification algorithms for use
% in paediatric pathological populations".

% Date: Modified for single file analysis without force plate
% Original: Sailee Sansgiri & Rosa Visscher

close all;
clc;

%% Initial set-up - MODIFY THESE PATHS
btkFolder     = 'D:\software\btk_0.3.0_Win7_MatlabR2009b_64bit\btk';%add location were btk folder from biomechanical toolkit is saved
addpath(btkFolder);

% 添加函数文件路径（假设所有f_*.m文件在当前文件夹）
addpath('D:\motiondata\HEAT_gait-event-detection\AGED_gait-event-detection-main');  % 当前文件夹，如果函数文件在别处请修改

% 设置C3D文件完整路径
c3dFilePath = 'D:\motiondata\50_StrokePiG\TVC03\BWA7.c3d';

%% 手动设置参数
side = 'Right';          % 'Left' 或 'Right'，根据实际情况设置

%% 检查文件是否存在
if ~exist(c3dFilePath, 'file')
    error('C3D文件不存在，请检查路径: %s', c3dFilePath);
end

%% 分析单个文件
fprintf('正在分析文件: %s\n', c3dFilePath);

% 读取C3D文件
btkData = btkReadAcquisition(c3dFilePath);
btkClearEvents(btkData);
metadata = btkGetMetaData(btkData);
ff = btkGetFirstFrame(btkData);
Markers = btkGetMarkers(btkData);
f = btkGetPointFrequency(btkData);
n = btkGetPointFrameNumber(btkData);

fprintf('采样频率: %.1f Hz\n', f);
fprintf('数据帧数: %d\n', n);

%% 设置标记点名称（基于您的可用标记点）
if strcmpi(side, 'Left')
    heelMarkerName = 'LHEE';
    toeMarkerName = 'LTOE';
    ankleMarkerName = 'LANK';  % 用于替代部分分析
else
    heelMarkerName = 'RHEE';
    toeMarkerName = 'RTOE';
    ankleMarkerName = 'RANK';  % 用于替代部分分析
end

% 使用骨盆标记点替代SACR
% 计算骨盆中心作为SACR的替代
LASI_data = Markers.LASI;
RASI_data = Markers.RASI;
LPSI_data = Markers.LPSI;
RPSI_data = Markers.RPSI;

% 计算骨盆中心（SACR替代）
SACR_substitute = (LASI_data + RASI_data + LPSI_data + RPSI_data) / 4;

fprintf('使用骨盆四个标记点的平均值替代SACR标记点\n');
fprintf('分析侧别: %s\n', side);
fprintf('使用标记点: %s (足跟), %s (足趾)\n', heelMarkerName, toeMarkerName);

%% 由于没有测力台数据，设置参考事件为NaN
mFS = NaN;
mFO = NaN;
mFS_ms = NaN;
mFO_ms = NaN;

fprintf('注意：没有使用测力台参考数据\n');

%% 步行方向矫正
% 使用骨盆中心数据进行方向判断
dir_i = abs(SACR_substitute(end, 1) - SACR_substitute(1, 1));
dir_j = abs(SACR_substitute(end, 2) - SACR_substitute(1, 2));

walkdir = 1;  % x is walkdir
if (dir_i < dir_j)
    walkdir = 2;  % y is walkdir
end

% 确定方向的正负
sgn = sign(SACR_substitute(end, walkdir) - SACR_substitute(1, walkdir));
walkdir = walkdir * sgn;

% 创建包含替代SACR的标记点结构
Markers_with_SACR = Markers;
Markers_with_SACR.SACR = SACR_substitute;

[Markers_Corrected] = f_rotCoordinateSystem(Markers_with_SACR, walkdir, 1);
gaitAxis = 3;
verticalAxis = 3;

%% 滤波标记点和预处理
[B, A] = butter(4, 6/(f/2), 'low');

% 滤波处理
filtheelmarker = filtfilt(B, A, Markers_Corrected.(heelMarkerName));
filttoemarker = filtfilt(B, A, Markers_Corrected.(toeMarkerName));
filtsacrmarker = filtfilt(B, A, Markers_Corrected.SACR);  % 使用替代的SACR
filtLPSI = filtfilt(B, A, Markers_Corrected.LPSI);
filtRPSI = filtfilt(B, A, Markers_Corrected.RPSI);
filtLASI = filtfilt(B, A, Markers_Corrected.LASI);
filtRASI = filtfilt(B, A, Markers_Corrected.RASI);

ysacr = filtsacrmarker(:, gaitAxis);
zsacr = filtsacrmarker(:, verticalAxis);
yheel = filtheelmarker(:, gaitAxis);
ytoe = filttoemarker(:, gaitAxis);

% 确定大致步行速度
[vel, time] = f_approxVelocity(ysacr, zsacr, f);
vel2 = vel/100;
fprintf('估计步行速度: %.2f m/s\n', vel2);

%% 运动学分析
% 计算标记点速度
Hvelocity_horizontal = zeros(n-1, 1);
Fvelocity_horizontal = zeros(n-1, 1);
Hvelocity_vertical = zeros(n-1, 1);
Fvelocity_vertical = zeros(n-1, 1);
Hvelocity_sagittal = zeros(n-1, 1);
Fvelocity_sagittal = zeros(n-1, 1);

for t = 1:n-1
    Hvelocity_sagittal(t) = sqrt((filtheelmarker(t+1,gaitAxis)- filtheelmarker(t,gaitAxis))^2+(filtheelmarker(t+1,verticalAxis)- filtheelmarker(t,verticalAxis))^2)/(1/f);
    Fvelocity_sagittal(t) = sqrt((filttoemarker(t+1,gaitAxis)- filttoemarker(t,gaitAxis)).^2+(filttoemarker(t+1,verticalAxis)- filttoemarker(t,verticalAxis)).^2)/(1/f);
    Hvelocity_horizontal(t) = (filtheelmarker(t+1,gaitAxis)-filtheelmarker(t,gaitAxis))/(1/f);
    Fvelocity_horizontal(t) = (filttoemarker(t+1,gaitAxis)-filttoemarker(t,gaitAxis))/(1/f);
    Hvelocity_vertical(t) = (filtheelmarker(t+1,verticalAxis)-filtheelmarker(t,verticalAxis))/(1/f);
    Fvelocity_vertical(t) = (filttoemarker(t+1,verticalAxis)-filttoemarker(t,verticalAxis))/(1/f);
end

%% 应用各种算法（无测力台参考）
fprintf('\n正在应用各种步态事件检测算法...\n');

% 1. Zeni算法
try
    [eFO_zeni_frame, eFS_zeni_frame, eFO_zeni_ms, eFS_zeni_ms] = f_zeni_event(f, yheel, ytoe, ysacr, mFO, mFS);
    fprintf('Zeni算法 - FS: 第%d帧 (%.1f ms), FO: 第%d帧 (%.1f ms)\n', eFS_zeni_frame, eFS_zeni_ms, eFO_zeni_frame, eFO_zeni_ms);
catch ME
    fprintf('Zeni算法执行失败: %s\n', ME.message);
    eFS_zeni_frame = NaN; eFO_zeni_frame = NaN;
    eFS_zeni_ms = NaN; eFO_zeni_ms = NaN;
end

% 2. Ghoussayni算法
try
    % 需要重新整形数据以匹配函数期望的格式
    heel_3d = reshape(filtheelmarker, [size(filtheelmarker,1), size(filtheelmarker,2), 1]);
    toe_3d = reshape(filttoemarker, [size(filttoemarker,1), size(filttoemarker,2), 1]);
    
    [FS_G, FO_G] = f_Ghoussayni_500(heel_3d, toe_3d, gaitAxis, verticalAxis, n, f);
    [eFO_G_frame, eFS_G_frame, eFO_G_ms, eFS_G_ms] = f_mG_event(FO_G, FS_G, mFO, mFS, f);
    fprintf('Ghoussayni算法 - FS: 第%d帧 (%.1f ms), FO: 第%d帧 (%.1f ms)\n', eFS_G_frame, eFS_G_ms, eFO_G_frame, eFO_G_ms);
catch ME
    fprintf('Ghoussayni算法执行失败: %s\n', ME.message);
    eFS_G_frame = NaN; eFO_G_frame = NaN;
    eFS_G_ms = NaN; eFO_G_ms = NaN;
end

% 3. 修改的Ghoussayni算法
try
    [FS_mG, FO_mG] = f_Ghoussayni_variablethreshold(heel_3d, toe_3d, gaitAxis, verticalAxis, n, f, vel2);
    [eFO_mG_frame, eFS_mG_frame, eFO_mG_ms, eFS_mG_ms] = f_mG_event(FO_mG, FS_mG, mFO, mFS, f);
    fprintf('修改Ghoussayni算法 - FS: 第%d帧 (%.1f ms), FO: 第%d帧 (%.1f ms)\n', eFS_mG_frame, eFS_mG_ms, eFO_mG_frame, eFO_mG_ms);
catch ME
    fprintf('修改Ghoussayni算法执行失败: %s\n', ME.message);
    eFS_mG_frame = NaN; eFO_mG_frame = NaN;
    eFS_mG_ms = NaN; eFO_mG_ms = NaN;
end

% 4. Desailly算法
try
    [B, A] = butter(4, (7/(f/2)));
    filttoemarker_d = filtfilt(B, A, Markers_Corrected.(toeMarkerName));
    filtheelmarker_d = filtfilt(B, A, Markers_Corrected.(heelMarkerName));
    fhm2 = filttoemarker_d(:, gaitAxis);
    fhm = filtheelmarker_d(:, gaitAxis);
    [z, p, k] = butter(4, 0.5/(f/2), 'high');
    [sos, g] = zp2sos(z, p, k);
    L_toe_high = filtfilt(sos, g, fhm2);
    L_heel_high = filtfilt(sos, g, fhm);
    
    [location_d_TO, index_d_TO] = findpeaks(-L_toe_high);
    [location_d_FS, index_d_FS] = findpeaks(L_heel_high);
    
    [eFS_D_frame, eFS_D_ms] = f_desailly_event(L_heel_high, mFS, f, location_d_FS, index_d_FS);
    [eFO_D_frame, eFO_D_ms] = f_desailly_event(L_heel_high, mFO, f, location_d_TO, index_d_TO);
    fprintf('Desailly算法 - FS: 第%d帧 (%.1f ms), FO: 第%d帧 (%.1f ms)\n', eFS_D_frame, eFS_D_ms, eFO_D_frame, eFO_D_ms);
catch ME
    fprintf('Desailly算法执行失败: %s\n', ME.message);
    eFS_D_frame = NaN; eFO_D_frame = NaN;
    eFS_D_ms = NaN; eFO_D_ms = NaN;
end

% 5. O'Connor算法
try
    zheel = filtheelmarker_d(:, verticalAxis);
    ztoe = filttoemarker_d(:, verticalAxis);
    footcentre = (zheel + ztoe) / 2;
    velfootcentre = diff(footcentre) * f;
    
    [eFO_oconnor_frame, eFS_oconnor_frame, eFO_oconnor_ms, eFS_oconnor_ms] = f_oConnor_event(f, zheel, velfootcentre, mFS, mFO);
    fprintf('O''Connor算法 - FS: 第%d帧 (%.1f ms), FO: 第%d帧 (%.1f ms)\n', eFS_oconnor_frame, eFS_oconnor_ms, eFO_oconnor_frame, eFO_oconnor_ms);
catch ME
    fprintf('O''Connor算法执行失败: %s\n', ME.message);
    eFS_oconnor_frame = NaN; eFO_oconnor_frame = NaN;
    eFS_oconnor_ms = NaN; eFO_oconnor_ms = NaN;
end

%% 简化的步态事件检测（基于速度阈值）
% 添加一个简单的速度阈值方法作为参考
fprintf('\n应用简化速度阈值方法...\n');
try
    % 计算足跟水平速度
    heel_horiz_vel = abs(Hvelocity_horizontal);
    toe_horiz_vel = abs(Fvelocity_horizontal);
    
    % 使用低速度阈值检测足跟着地
    vel_threshold_fs = 100; % mm/s
    vel_threshold_fo = 200; % mm/s
    
    % 找到速度低于阈值的点（足跟着地）
    fs_candidates = find(heel_horiz_vel < vel_threshold_fs);
    % 找到速度高于阈值的点（足趾离地）
    fo_candidates = find(toe_horiz_vel > vel_threshold_fo);
    
    if ~isempty(fs_candidates)
        % 找到第一个和最后一个低速度点
        simple_FS1 = fs_candidates(1);
        if length(fs_candidates) > 50
            simple_FS2 = fs_candidates(end);
        else
            simple_FS2 = NaN;
        end
    else
        simple_FS1 = NaN;
        simple_FS2 = NaN;
    end
    
    if ~isempty(fo_candidates)
        simple_FO1 = fo_candidates(1);
        if length(fo_candidates) > 50
            simple_FO2 = fo_candidates(end);
        else
            simple_FO2 = NaN;
        end
    else
        simple_FO1 = NaN;
        simple_FO2 = NaN;
    end
    
    fprintf('简化方法 - 第1个FS: 第%d帧, 第1个FO: 第%d帧\n', simple_FS1, simple_FO1);
    if ~isnan(simple_FS2)
        fprintf('简化方法 - 第2个FS: 第%d帧, 第2个FO: 第%d帧\n', simple_FS2, simple_FO2);
    end
catch ME
    fprintf('简化方法执行失败: %s\n', ME.message);
    simple_FS1 = NaN; simple_FO1 = NaN;
end

%% 结果汇总和算法间比较分析
fprintf('\n=== 步态事件检测结果汇总 ===\n');
fprintf('文件: %s\n', c3dFilePath);
fprintf('侧别: %s\n', side);
fprintf('采样频率: %.1f Hz\n', f);
fprintf('步行速度: %.2f m/s\n', vel2);
fprintf('\n主要算法检测结果:\n');
fprintf('%-20s - FS: %4d帧 (%6.1f ms), FO: %4d帧 (%6.1f ms)\n', 'Zeni算法', eFS_zeni_frame, eFS_zeni_ms, eFO_zeni_frame, eFO_zeni_ms);
fprintf('%-20s - FS: %4d帧 (%6.1f ms), FO: %4d帧 (%6.1f ms)\n', 'Ghoussayni算法', eFS_G_frame, eFS_G_ms, eFO_G_frame, eFO_G_ms);
fprintf('%-20s - FS: %4d帧 (%6.1f ms), FO: %4d帧 (%6.1f ms)\n', '修改Ghoussayni', eFS_mG_frame, eFS_mG_ms, eFO_mG_frame, eFO_mG_ms);
fprintf('%-20s - FS: %4d帧 (%6.1f ms), FO: %4d帧 (%6.1f ms)\n', 'Desailly算法', eFS_D_frame, eFS_D_ms, eFO_D_frame, eFO_D_ms);
fprintf('%-20s - FS: %4d帧 (%6.1f ms), FO: %4d帧 (%6.1f ms)\n', 'O''Connor算法', eFS_oconnor_frame, eFS_oconnor_ms, eFO_oconnor_frame, eFO_oconnor_ms);

%% 算法间比较分析
fprintf('\n=== 算法间比较分析 ===\n');

% 收集有效结果
valid_fs_results = [];
valid_fo_results = [];
valid_methods_fs = {};
valid_methods_fo = {};

fs_all = [eFS_zeni_frame, eFS_G_frame, eFS_mG_frame, eFS_D_frame, eFS_oconnor_frame];
fo_all = [eFO_zeni_frame, eFO_G_frame, eFO_mG_frame, eFO_D_frame, eFO_oconnor_frame];
method_names = {'Zeni', 'Ghoussayni', 'Modified_Ghoussayni', 'Desailly', 'OConnor'};

for i = 1:length(fs_all)
    if ~isnan(fs_all(i))
        valid_fs_results(end+1) = fs_all(i);
        valid_methods_fs{end+1} = method_names{i};
    end
    if ~isnan(fo_all(i))
        valid_fo_results(end+1) = fo_all(i);
        valid_methods_fo{end+1} = method_names{i};
    end
end

% 足跟着地事件比较
if length(valid_fs_results) > 1
    fprintf('足跟着地(FS)事件比较:\n');
    fs_mean = mean(valid_fs_results);
    fs_std = std(valid_fs_results);
    fs_range = max(valid_fs_results) - min(valid_fs_results);
    
    fprintf('  有效检测算法数: %d/%d\n', length(valid_fs_results), length(method_names));
    fprintf('  平均值: %.1f帧 (%.1f ms)\n', fs_mean, fs_mean/f*1000);
    fprintf('  标准差: %.1f帧 (%.1f ms)\n', fs_std, fs_std/f*1000);
    fprintf('  范围: %.1f帧 (%.1f ms)\n', fs_range, fs_range/f*1000);
    
    fprintf('  各算法与平均值的差异:\n');
    for i = 1:length(valid_fs_results)
        diff_frames = valid_fs_results(i) - fs_mean;
        diff_ms = diff_frames / f * 1000;
        fprintf('    %s: %+.1f帧 (%+.1f ms)\n', valid_methods_fs{i}, diff_frames, diff_ms);
    end
else
    fprintf('足跟着地(FS): 只有%d个算法成功检测\n', length(valid_fs_results));
end

% 足趾离地事件比较
if length(valid_fo_results) > 1
    fprintf('\n足趾离地(FO)事件比较:\n');
    fo_mean = mean(valid_fo_results);
    fo_std = std(valid_fo_results);
    fo_range = max(valid_fo_results) - min(valid_fo_results);
    
    fprintf('  有效检测算法数: %d/%d\n', length(valid_fo_results), length(method_names));
    fprintf('  平均值: %.1f帧 (%.1f ms)\n', fo_mean, fo_mean/f*1000);
    fprintf('  标准差: %.1f帧 (%.1f ms)\n', fo_std, fo_std/f*1000);
    fprintf('  范围: %.1f帧 (%.1f ms)\n', fo_range, fo_range/f*1000);
    
    fprintf('  各算法与平均值的差异:\n');
    for i = 1:length(valid_fo_results)
        diff_frames = valid_fo_results(i) - fo_mean;
        diff_ms = diff_frames / f * 1000;
        fprintf('    %s: %+.1f帧 (%+.1f ms)\n', valid_methods_fo{i}, diff_frames, diff_ms);
    end
else
    fprintf('足趾离地(FO): 只有%d个算法成功检测\n', length(valid_fo_results));
end

% 计算步态相关参数
if length(valid_fs_results) > 0 && length(valid_fo_results) > 0
    fprintf('\n步态时间参数:\n');
    % 使用最接近平均值的算法结果
    if length(valid_fs_results) > 1
        [~, best_fs_idx] = min(abs(valid_fs_results - mean(valid_fs_results)));
        best_fs = valid_fs_results(best_fs_idx);
    else
        best_fs = valid_fs_results(1);
    end
    
    if length(valid_fo_results) > 1
        [~, best_fo_idx] = min(abs(valid_fo_results - mean(valid_fo_results)));
        best_fo = valid_fo_results(best_fo_idx);
    else
        best_fo = valid_fo_results(1);
    end
    
    if best_fo > best_fs
        stance_time = (best_fo - best_fs) / f * 1000;
        fprintf('  支撑相时间: %.1f ms (第%d帧到第%d帧)\n', stance_time, best_fs, best_fo);
    end
end

%% 创建结果表格
results = table();
results.Method = {'Zeni'; 'Ghoussayni'; 'Modified_Ghoussayni'; 'Desailly'; 'OConnor'};
results.FS_Frame = [eFS_zeni_frame; eFS_G_frame; eFS_mG_frame; eFS_D_frame; eFS_oconnor_frame];
results.FO_Frame = [eFO_zeni_frame; eFO_G_frame; eFO_mG_frame; eFO_D_frame; eFO_oconnor_frame];
results.FS_ms = [eFS_zeni_ms; eFS_G_ms; eFS_mG_ms; eFS_D_ms; eFS_oconnor_ms];
results.FO_ms = [eFO_zeni_ms; eFO_G_ms; eFO_mG_ms; eFO_D_ms; eFO_oconnor_ms];

fprintf('\n详细结果表格:\n');
disp(results);

%% 创建算法一致性分析图
if length(valid_fs_results) > 1 || length(valid_fo_results) > 1
    figure('Position', [150, 150, 800, 400]);
    
    if length(valid_fs_results) > 1 && length(valid_fo_results) > 1
        subplot(1,2,1);
        scatter(1:length(valid_fs_results), valid_fs_results, 100, 'filled', 'blue'); hold on;
        plot([1, length(valid_fs_results)], [mean(valid_fs_results), mean(valid_fs_results)], 'r--', 'LineWidth', 2);
        plot([1, length(valid_fs_results)], [mean(valid_fs_results)+std(valid_fs_results), mean(valid_fs_results)+std(valid_fs_results)], 'r:', 'LineWidth', 1);
        plot([1, length(valid_fs_results)], [mean(valid_fs_results)-std(valid_fs_results), mean(valid_fs_results)-std(valid_fs_results)], 'r:', 'LineWidth', 1);
        
        % 添加数值标签
        for i = 1:length(valid_fs_results)
            text(i, valid_fs_results(i)+5, sprintf('%.0f', valid_fs_results(i)), ...
                'HorizontalAlignment', 'center', 'FontSize', 10);
        end
        
        set(gca, 'XTick', 1:length(valid_methods_fs), 'XTickLabel', valid_methods_fs);
        xlabel('算法'); ylabel('帧数'); 
        title('足跟着地(FS)检测结果一致性');
        legend('检测结果', '平均值', '±1标准差', 'Location', 'best');
        grid on;
        
        subplot(1,2,2);
        scatter(1:length(valid_fo_results), valid_fo_results, 100, 'filled', 'red'); hold on;
        plot([1, length(valid_fo_results)], [mean(valid_fo_results), mean(valid_fo_results)], 'b--', 'LineWidth', 2);
        plot([1, length(valid_fo_results)], [mean(valid_fo_results)+std(valid_fo_results), mean(valid_fo_results)+std(valid_fo_results)], 'b:', 'LineWidth', 1);
        plot([1, length(valid_fo_results)], [mean(valid_fo_results)-std(valid_fo_results), mean(valid_fo_results)-std(valid_fo_results)], 'b:', 'LineWidth', 1);
        
        % 添加数值标签
        for i = 1:length(valid_fo_results)
            text(i, valid_fo_results(i)+5, sprintf('%.0f', valid_fo_results(i)), ...
                'HorizontalAlignment', 'center', 'FontSize', 10);
        end
        
        set(gca, 'XTick', 1:length(valid_methods_fo), 'XTickLabel', valid_methods_fo);
        xlabel('算法'); ylabel('帧数'); 
        title('足趾离地(FO)检测结果一致性');
        legend('检测结果', '平均值', '±1标准差', 'Location', 'best');
        grid on;
    else
        % 如果只有一种事件有多个检测结果
        if length(valid_fs_results) > 1
            scatter(1:length(valid_fs_results), valid_fs_results, 100, 'filled', 'blue'); hold on;
            plot([1, length(valid_fs_results)], [mean(valid_fs_results), mean(valid_fs_results)], 'r--', 'LineWidth', 2);
            set(gca, 'XTick', 1:length(valid_methods_fs), 'XTickLabel', valid_methods_fs);
            title('足跟着地(FS)检测结果一致性');
        else
            scatter(1:length(valid_fo_results), valid_fo_results, 100, 'filled', 'red'); hold on;
            plot([1, length(valid_fo_results)], [mean(valid_fo_results), mean(valid_fo_results)], 'b--', 'LineWidth', 2);
            set(gca, 'XTick', 1:length(valid_methods_fo), 'XTickLabel', valid_methods_fo);
            title('足趾离地(FO)检测结果一致性');
        end
        xlabel('算法'); ylabel('帧数'); 
        legend('检测结果', '平均值', 'Location', 'best');
        grid on;
    end
end

%% 绘制简单的可视化图形
figure('Position', [100, 100, 1200, 600]);

% 绘制足跟和足趾的水平位置
subplot(2,2,1);
plot(yheel, 'b-', 'LineWidth', 1.5); hold on;
plot(ytoe, 'r-', 'LineWidth', 1.5);
if ~isnan(eFS_zeni_frame)
    plot(eFS_zeni_frame, yheel(eFS_zeni_frame), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
end
if ~isnan(eFO_zeni_frame)
    plot(eFO_zeni_frame, ytoe(eFO_zeni_frame), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
end
xlabel('帧数'); ylabel('位置 (mm)'); 
title('足跟和足趾水平位置 (Zeni算法检测点)');
legend('足跟', '足趾', 'FS检测', 'FO检测');
grid on;

% 绘制速度曲线
subplot(2,2,2);
plot(abs(Hvelocity_horizontal), 'b-', 'LineWidth', 1); hold on;
plot(abs(Fvelocity_horizontal), 'r-', 'LineWidth', 1);
xlabel('帧数'); ylabel('速度 (mm/s)'); 
title('足跟和足趾水平速度');
legend('足跟速度', '足趾速度');
grid on;

% 绘制Zeni算法的差值信号
subplot(2,2,3);
tFS = yheel - ysacr;
tFO = ytoe - ysacr;
plot(tFS, 'b-', 'LineWidth', 1); hold on;
plot(tFO, 'r-', 'LineWidth', 1);
if ~isnan(eFS_zeni_frame)
    plot(eFS_zeni_frame, tFS(eFS_zeni_frame), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
end
if ~isnan(eFO_zeni_frame)
    plot(eFO_zeni_frame, tFO(eFO_zeni_frame), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
end
xlabel('帧数'); ylabel('相对位置 (mm)'); 
title('Zeni算法: 足部相对骨盆位置');
legend('足跟-骨盆', '足趾-骨盆', 'FS检测', 'FO检测');
grid on;

% 绘制检测结果比较
subplot(2,2,4);
methods = {'Zeni', 'Ghouss', 'mGhouss', 'Desailly', 'O''Connor'};
fs_results = [eFS_zeni_frame, eFS_G_frame, eFS_mG_frame, eFS_D_frame, eFS_oconnor_frame];
fo_results = [eFO_zeni_frame, eFO_G_frame, eFO_mG_frame, eFO_D_frame, eFO_oconnor_frame];

% 移除NaN值进行绘图
valid_fs = ~isnan(fs_results);
valid_fo = ~isnan(fo_results);

if any(valid_fs)
    bar(find(valid_fs), fs_results(valid_fs), 'b', 'FaceAlpha', 0.7); hold on;
end
if any(valid_fo)
    bar(find(valid_fo), fo_results(valid_fo), 'r', 'FaceAlpha', 0.7);
end

set(gca, 'XTick', 1:length(methods), 'XTickLabel', methods);
xlabel('算法'); ylabel('帧数'); 
title('各算法检测结果比较');
legend('足跟着地 (FS)', '足趾离地 (FO)');
grid on;

fprintf('\n分析完成！检查图形窗口查看可视化结果。\n');