%% 修复版步态事件检测 - 完整解决方案
% 包含修复的算法函数和改进的检测方法

close all;
clc;

%% Initial set-up
btkFolder = 'D:\VR_projects\motionanalysis\btk-0.3.0_Win7_MatlabR2009b_64bit';  % 修改为您的BTK安装路径
addpath(btkFolder);
addpath('D:\VR_projects\motionanalysis\AGED_gait-event-detection-main');

c3dFilePath = 'D:\VR_projects\motionanalysis\50_StrokePiG\TVC03\BWA7.c3d';
side = 'Right';

%% 读取数据
fprintf('=== 读取C3D文件 ===\n');
btkData = btkReadAcquisition(c3dFilePath);
btkClearEvents(btkData);
Markers = btkGetMarkers(btkData);
f = btkGetPointFrequency(btkData);
n = btkGetPointFrameNumber(btkData);

fprintf('文件: %s\n', c3dFilePath);
fprintf('采样频率: %.1f Hz\n', f);
fprintf('数据帧数: %d (%.2f 秒)\n', n, n/f);

%% 设置标记点
if strcmpi(side, 'Left')
    heelMarkerName = 'LHEE';
    toeMarkerName = 'LTOE';
else
    heelMarkerName = 'RHEE';
    toeMarkerName = 'RTOE';
end

%% 创建骨盆中心（替代SACR）
LASI_data = Markers.LASI;
RASI_data = Markers.RASI;
LPSI_data = Markers.LPSI;
RPSI_data = Markers.RPSI;
SACR_substitute = (LASI_data + RASI_data + LPSI_data + RPSI_data) / 4;

%% 坐标系矫正
displacement = SACR_substitute(end, :) - SACR_substitute(1, :);
dir_i = abs(displacement(1));
dir_j = abs(displacement(2));

walkdir = 1;
if (dir_i < dir_j)
    walkdir = 2;
end

sgn = sign(displacement(walkdir));
walkdir = walkdir * sgn;

% 简化的坐标矫正（如果原函数有问题）
Markers_Corrected = Markers;
Markers_Corrected.SACR = SACR_substitute;

% 如果有坐标矫正函数且工作正常
try
    [Markers_Corrected] = f_rotCoordinateSystem(Markers_Corrected, walkdir, 1);
    fprintf('坐标系矫正成功\n');
catch
    fprintf('使用原始坐标系\n');
    Markers_Corrected.SACR = SACR_substitute;
end

gaitAxis = 1;
verticalAxis = 3;

%% 数据滤波
[B, A] = butter(4, 6/(f/2), 'low');
filtheelmarker = filtfilt(B, A, Markers_Corrected.(heelMarkerName));
filttoemarker = filtfilt(B, A, Markers_Corrected.(toeMarkerName));
filtsacrmarker = filtfilt(B, A, Markers_Corrected.SACR);

ysacr = filtsacrmarker(:, gaitAxis);
zsacr = filtsacrmarker(:, verticalAxis);
yheel = filtheelmarker(:, gaitAxis);
ytoe = filttoemarker(:, gaitAxis);
zheel = filtheelmarker(:, verticalAxis);
ztoe = filttoemarker(:, verticalAxis);

%% checkpoint 1
% figure
% subplot(1,2,1)
% hh = Markers_Corrected.(heelMarkerName);
% plot(filtheelmarker(:,3));
% 
% subplot(1,2,2)
% plot(hh(:,3));

%% 重新计算步行速度
total_disp = sqrt(sum(displacement(1:2).^2));  % 水平位移
walking_time = n/f;
corrected_speed = total_disp / walking_time / 1000;  % m/s

fprintf('总水平位移: %.1f mm\n', total_disp);
fprintf('步行时间: %.2f 秒\n', walking_time);
fprintf('修正步行速度: %.2f m/s\n', corrected_speed);

% 使用合理的速度值
if corrected_speed > 0.3 && corrected_speed < 3.0
    vel2 = corrected_speed;
else
    vel2 = 1.2;  % 默认步行速度
    fprintf('使用默认步行速度: %.2f m/s\n', vel2);
end

%% 内嵌修复的Zeni算法
fprintf('\n=== 应用Zeni算法 ===\n');
try
    tFS = yheel - ysacr;
    tFO = ytoe - ysacr;
    
    [vFS_zeni, FSindex_zeni] = findpeaks(tFS, 'MinPeakDistance', f*0.3);  % 最小间隔0.3秒
    [vFO_zeni, FOindex_zeni] = findpeaks(-tFO, 'MinPeakDistance', f*0.3);
%% checkpoint 2
    % %% 可视化FS检测结果
    % fig_fs = plotZeniEvents(tFS, FSindex_zeni, vFS_zeni, ...
    %                        'EventType', 'FS', ...
    %                        'SamplingFreq', f, ...
    %                        'ShowTime', true);
    % fo_original_values = tFO(FOindex_zeni);            % ✅ 获取原始tFO的值（负数）
    % fig_fo = plotZeniEvents(tFO, FOindex_zeni, fo_original_values, ...
    %                    'EventType', 'FO', ...
    %                    'SamplingFreq', f, ...
    %                    'ShowTime', true);

%%
    fprintf('Zeni算法找到 %d 个FS峰值, %d 个FO峰值\n', length(FSindex_zeni), length(FOindex_zeni));
    
    if ~isempty(FSindex_zeni)
        eFS_zeni_frame = FSindex_zeni(1);
        eFS_zeni_ms = eFS_zeni_frame / f * 1000;
    else
        eFS_zeni_frame = NaN;
        eFS_zeni_ms = NaN;
    end
    
    if ~isempty(FOindex_zeni)
        eFO_zeni_frame = FOindex_zeni(1);
        eFO_zeni_ms = eFO_zeni_frame / f * 1000;
    else
        eFO_zeni_frame = NaN;
        eFO_zeni_ms = NaN;
    end
    
    fprintf('Zeni结果: FS=%d帧(%.1fms), FO=%d帧(%.1fms)\n', ...
        eFS_zeni_frame, eFS_zeni_ms, eFO_zeni_frame, eFO_zeni_ms);
    
catch ME
    fprintf('Zeni算法失败: %s\n', ME.message);
    eFS_zeni_frame = NaN; eFO_zeni_frame = NaN;
    eFS_zeni_ms = NaN; eFO_zeni_ms = NaN;
end

%% 内嵌修复的速度阈值算法
fprintf('\n=== 应用速度阈值算法 ===\n');
try
    % 计算2D速度（水平+垂直）
    heel_vel_2d = zeros(n-1, 1);
    toe_vel_2d = zeros(n-1, 1);
    
    for t = 1:n-1
        heel_vel_2d(t) = sqrt((filtheelmarker(t+1,gaitAxis) - filtheelmarker(t,gaitAxis))^2 + ...
                             (filtheelmarker(t+1,verticalAxis) - filtheelmarker(t,verticalAxis))^2) * f;
        toe_vel_2d(t) = sqrt((filttoemarker(t+1,gaitAxis) - filttoemarker(t,gaitAxis))^2 + ...
                            (filttoemarker(t+1,verticalAxis) - filttoemarker(t,verticalAxis))^2) * f;
    end
    
    % 自适应阈值基于速度分布
    heel_thresh = prctile(heel_vel_2d, 15);  % 足跟低速阈值
    toe_thresh = prctile(toe_vel_2d, 85);    % 足趾高速阈值
    
    fprintf('足跟速度范围: [%.1f, %.1f] mm/s, 阈值: %.1f mm/s\n', ...
        min(heel_vel_2d), max(heel_vel_2d), heel_thresh);
    fprintf('足趾速度范围: [%.1f, %.1f] mm/s, 阈值: %.1f mm/s\n', ...
        min(toe_vel_2d), max(toe_vel_2d), toe_thresh);
    
    % 寻找连续低速区域作为足跟着地
    fs_candidates = find(heel_vel_2d < heel_thresh);
    fo_candidates = find(toe_vel_2d > toe_thresh);
    
    % 分组连续的候选点
    if ~isempty(fs_candidates)
        fs_groups = [];
        current_group = fs_candidates(1);
        for i = 2:length(fs_candidates)
            if fs_candidates(i) - fs_candidates(i-1) > 5  % 间隔大于5帧认为是新组
                fs_groups(end+1) = round(mean(current_group));
                current_group = fs_candidates(i);
            else
                current_group(end+1) = fs_candidates(i);
            end
        end
        fs_groups(end+1) = round(mean(current_group));
        
        eFS_vel_frame = fs_groups(1);
        eFS_vel_ms = eFS_vel_frame / f * 1000;
    else
        eFS_vel_frame = NaN;
        eFS_vel_ms = NaN;
    end
    
    % 分组连续的FO候选点
    if ~isempty(fo_candidates)
        fo_groups = [];
        current_group = fo_candidates(1);
        for i = 2:length(fo_candidates)
            if fo_candidates(i) - fo_candidates(i-1) > 5
                fo_groups(end+1) = round(mean(current_group));
                current_group = fo_candidates(i);
            else
                current_group(end+1) = fo_candidates(i);
            end
        end
        fo_groups(end+1) = round(mean(current_group));
        
        eFO_vel_frame = fo_groups(1);
        eFO_vel_ms = eFO_vel_frame / f * 1000;
    else
        eFO_vel_frame = NaN;
        eFO_vel_ms = NaN;
    end
    
    fprintf('速度阈值结果: FS=%d帧(%.1fms), FO=%d帧(%.1fms)\n', ...
        eFS_vel_frame, eFS_vel_ms, eFO_vel_frame, eFO_vel_ms);
    
catch ME
    fprintf('速度阈值算法失败: %s\n', ME.message);
    eFS_vel_frame = NaN; eFO_vel_frame = NaN;
    eFS_vel_ms = NaN; eFO_vel_ms = NaN;
end

%% 内嵌O'Connor风格算法
fprintf('\n=== 应用O''Connor风格算法 ===\n');
try
    % 计算足部中心的垂直速度
    foot_center_z = (zheel + ztoe) / 2;
    foot_center_vel_z = diff(foot_center_z) * f;
    
    % 寻找垂直速度的峰值和谷值
    [~, fs_oconnor_idx] = findpeaks(-foot_center_vel_z, 'MinPeakDistance', f*0.3);
    [~, fo_oconnor_idx] = findpeaks(foot_center_vel_z, 'MinPeakDistance', f*0.3);
    
    % 添加足跟高度约束
    heel_height_thresh = min(zheel) + 0.3 * (max(zheel) - min(zheel));
    
    % 过滤FS候选点（足跟不能太高）
    valid_fs = [];
    for i = 1:length(fs_oconnor_idx)
        if fs_oconnor_idx(i) <= length(zheel) && zheel(fs_oconnor_idx(i)) < heel_height_thresh
            valid_fs(end+1) = fs_oconnor_idx(i);
        end
    end
    
    if ~isempty(valid_fs)
        eFS_oconnor_frame = valid_fs(1);
        eFS_oconnor_ms = eFS_oconnor_frame / f * 1000;
    else
        eFS_oconnor_frame = NaN;
        eFS_oconnor_ms = NaN;
    end
    
    if ~isempty(fo_oconnor_idx)
        eFO_oconnor_frame = fo_oconnor_idx(1);
        eFO_oconnor_ms = eFO_oconnor_frame / f * 1000;
    else
        eFO_oconnor_frame = NaN;
        eFO_oconnor_ms = NaN;
    end
    
    fprintf('O''Connor结果: FS=%d帧(%.1fms), FO=%d帧(%.1fms)\n', ...
        eFS_oconnor_frame, eFS_oconnor_ms, eFO_oconnor_frame, eFO_oconnor_ms);
    
catch ME
    fprintf('O''Connor算法失败: %s\n', ME.message);
    eFS_oconnor_frame = NaN; eFO_oconnor_frame = NaN;
    eFS_oconnor_ms = NaN; eFO_oconnor_ms = NaN;
end

%% 简化的Desailly风格算法
fprintf('\n=== 应用简化Desailly算法 ===\n');
try
    % 高通滤波
    [z, p, k] = butter(4, 1/(f/2), 'high');  % 1Hz高通滤波
    [sos, g] = zp2sos(z, p, k);
    
    heel_highpass = filtfilt(sos, g, yheel);
    toe_highpass = filtfilt(sos, g, ytoe);
    
    % 寻找峰值
    [~, fs_desailly_idx] = findpeaks(heel_highpass, 'MinPeakDistance', f*0.3);
    [~, fo_desailly_idx] = findpeaks(-toe_highpass, 'MinPeakDistance', f*0.3);
    
    if ~isempty(fs_desailly_idx)
        eFS_desailly_frame = fs_desailly_idx(1);
        eFS_desailly_ms = eFS_desailly_frame / f * 1000;
    else
        eFS_desailly_frame = NaN;
        eFS_desailly_ms = NaN;
    end
    
    if ~isempty(fo_desailly_idx)
        eFO_desailly_frame = fo_desailly_idx(1);
        eFO_desailly_ms = eFO_desailly_frame / f * 1000;
    else
        eFO_desailly_frame = NaN;
        eFO_desailly_ms = NaN;
    end
    
    fprintf('Desailly结果: FS=%d帧(%.1fms), FO=%d帧(%.1fms)\n', ...
        eFS_desailly_frame, eFS_desailly_ms, eFO_desailly_frame, eFO_desailly_ms);
    
catch ME
    fprintf('Desailly算法失败: %s\n', ME.message);
    eFS_desailly_frame = NaN; eFO_desailly_frame = NaN;
    eFS_desailly_ms = NaN; eFO_desailly_ms = NaN;
end

%% 结果汇总
fprintf('\n=== 步态事件检测结果汇总 ===\n');
fprintf('文件: %s\n', c3dFilePath);
fprintf('侧别: %s\n', side);
fprintf('采样频率: %.1f Hz\n', f);
fprintf('修正步行速度: %.2f m/s\n', vel2);
fprintf('\n各算法检测结果:\n');
fprintf('%-15s - FS: %4d帧 (%6.1f ms), FO: %4d帧 (%6.1f ms)\n', 'Zeni', eFS_zeni_frame, eFS_zeni_ms, eFO_zeni_frame, eFO_zeni_ms);
fprintf('%-15s - FS: %4d帧 (%6.1f ms), FO: %4d帧 (%6.1f ms)\n', '速度阈值', eFS_vel_frame, eFS_vel_ms, eFO_vel_frame, eFO_vel_ms);
fprintf('%-15s - FS: %4d帧 (%6.1f ms), FO: %4d帧 (%6.1f ms)\n', 'O''Connor', eFS_oconnor_frame, eFS_oconnor_ms, eFO_oconnor_frame, eFO_oconnor_ms);
fprintf('%-15s - FS: %4d帧 (%6.1f ms), FO: %4d帧 (%6.1f ms)\n', 'Desailly', eFS_desailly_frame, eFS_desailly_ms, eFO_desailly_frame, eFO_desailly_ms);

%% 算法一致性分析
fs_results = [eFS_zeni_frame, eFS_vel_frame, eFS_oconnor_frame, eFS_desailly_frame];
fo_results = [eFO_zeni_frame, eFO_vel_frame, eFO_oconnor_frame, eFO_desailly_frame];
method_names = {'Zeni', 'Velocity', 'O''Connor', 'Desailly'};

valid_fs = fs_results(~isnan(fs_results));
valid_fo = fo_results(~isnan(fo_results));

fprintf('\n=== 算法一致性分析 ===\n');
if length(valid_fs) > 1
    fs_mean = mean(valid_fs);
    fs_std = std(valid_fs);
    fprintf('FS检测 - 成功算法: %d/%d, 平均: %.1f帧, 标准差: %.1f帧\n', ...
        length(valid_fs), length(fs_results), fs_mean, fs_std);
else
    fprintf('FS检测 - 只有%d个算法成功\n', length(valid_fs));
end

if length(valid_fo) > 1
    fo_mean = mean(valid_fo);
    fo_std = std(valid_fo);
    fprintf('FO检测 - 成功算法: %d/%d, 平均: %.1f帧, 标准差: %.1f帧\n', ...
        length(valid_fo), length(fo_results), fo_mean, fo_std);
else
    fprintf('FO检测 - 只有%d个算法成功\n', length(valid_fo));
end

%% 可视化结果
figure('Position', [100, 100, 1400, 900]);

% 子图1: 原始轨迹和Zeni信号
subplot(2,3,1);
zheel = filtheelmarker(:, verticalAxis);  % 足跟垂直位置
ztoe = filttoemarker(:, verticalAxis);    % 足趾垂直位置
zsacr = filtsacrmarker(:, verticalAxis);  % 骨盆垂直位置

plot(zheel, 'b-', 'LineWidth', 1.5); hold on;
plot(ztoe, 'r-', 'LineWidth', 1.5);
plot(zsacr, 'g--', 'LineWidth', 1);

% 添加检测点（如果需要在垂直位置图上显示）
if ~isnan(eFS_zeni_frame)
    plot(eFS_zeni_frame, zheel(eFS_zeni_frame), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
end
if ~isnan(eFO_zeni_frame)
    plot(eFO_zeni_frame, ztoe(eFO_zeni_frame), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
end

xlabel('帧数'); 
ylabel('Z位置 (mm)');  % 修改Y轴标签
title('垂直位置轨迹 + Zeni检测');  % 修改标题
legend('足跟', '足趾', '骨盆中心', 'FS检测', 'FO检测');
grid on;

% 子图2: Zeni算法信号
subplot(2,3,2);
if exist('tFS', 'var')
    plot(tFS, 'b-', 'LineWidth', 1.5); hold on;
    plot(-tFO, 'r-', 'LineWidth', 1.5);
    if ~isnan(eFS_zeni_frame)
        plot(eFS_zeni_frame, tFS(eFS_zeni_frame), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
    end
    if ~isnan(eFO_zeni_frame)
        plot(eFO_zeni_frame, -tFO(eFO_zeni_frame), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    end
    xlabel('帧数'); ylabel('相对位置 (mm)');
    title('Zeni算法信号');
    legend('足跟-骨盆', '-(足趾-骨盆)', 'FS峰值', 'FO峰值');
    grid on;
end

% 子图3: 速度信号
subplot(2,3,3);
if exist('heel_vel_2d', 'var')
    plot(heel_vel_2d, 'b-', 'LineWidth', 1); hold on;
    plot(toe_vel_2d, 'r-', 'LineWidth', 1);
    plot([1, length(heel_vel_2d)], [heel_thresh, heel_thresh], 'b--', 'LineWidth', 2);
    plot([1, length(toe_vel_2d)], [toe_thresh, toe_thresh], 'r--', 'LineWidth', 2);
    if ~isnan(eFS_vel_frame)
        plot(eFS_vel_frame, heel_vel_2d(eFS_vel_frame), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
    end
    if ~isnan(eFO_vel_frame)
        plot(eFO_vel_frame, toe_vel_2d(eFO_vel_frame), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    end
    xlabel('帧数'); ylabel('速度 (mm/s)');
    title('2D速度 + 阈值检测');
    legend('足跟速度', '足趾速度', 'FS阈值', 'FO阈值', 'FS检测', 'FO检测');
    grid on;
end

% 子图4: 垂直位置
subplot(2,3,4);
plot(zheel, 'b-', 'LineWidth', 1.5); hold on;
plot(ztoe, 'r-', 'LineWidth', 1.5);
if exist('foot_center_z', 'var')
    plot(foot_center_z, 'g-', 'LineWidth', 1);
end
if ~isnan(eFS_oconnor_frame)
    plot(eFS_oconnor_frame, zheel(eFS_oconnor_frame), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
end
if ~isnan(eFO_oconnor_frame)
    plot(eFO_oconnor_frame, ztoe(eFO_oconnor_frame), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
end
xlabel('帧数'); ylabel('Z位置 (mm)');
title('垂直位置 + O''Connor检测');
legend('足跟', '足趾', '足部中心', 'FS', 'FO');
grid on;

% 子图5: 所有算法比较
subplot(2,3,5);
x_pos = 1:length(method_names);
bar_width = 0.35;

% 处理NaN值用于显示
fs_display = fs_results;
fo_display = fo_results;
fs_display(isnan(fs_results)) = 0;
fo_display(isnan(fo_results)) = 0;

h1 = bar(x_pos - bar_width/2, fs_display, bar_width, 'FaceColor', [0.3 0.7 1]); hold on;
h2 = bar(x_pos + bar_width/2, fo_display, bar_width, 'FaceColor', [1 0.5 0.3]);

% 标记有效数值
for i = 1:length(fs_results)
    if ~isnan(fs_results(i))
        text(i - bar_width/2, fs_results(i) + 20, num2str(fs_results(i)), ...
             'HorizontalAlignment', 'center', 'FontSize', 8);
    else
        text(i - bar_width/2, 50, 'NaN', 'HorizontalAlignment', 'center', 'FontSize', 8, 'Color', 'red');
    end
    
    if ~isnan(fo_results(i))
        text(i + bar_width/2, fo_results(i) + 20, num2str(fo_results(i)), ...
             'HorizontalAlignment', 'center', 'FontSize', 8);
    else
        text(i + bar_width/2, 50, 'NaN', 'HorizontalAlignment', 'center', 'FontSize', 8, 'Color', 'red');
    end
end

set(gca, 'XTick', x_pos, 'XTickLabel', method_names);
xlabel('算法'); ylabel('帧数');
title('各算法检测结果比较');
legend([h1, h2], {'足跟着地 (FS)', '足趾离地 (FO)'});
grid on;
max_val = max([valid_fs, valid_fo, 100]);
ylim([0, max_val + 100]);

% 子图6: 结果统计
subplot(2,3,6);
text(0.1, 0.9, sprintf('成功检测算法数量:'), 'FontSize', 12, 'FontWeight', 'bold');
text(0.1, 0.8, sprintf('  FS: %d/%d', length(valid_fs), length(fs_results)), 'FontSize', 11);
text(0.1, 0.7, sprintf('  FO: %d/%d', length(valid_fo), length(fo_results)), 'FontSize', 11);

if length(valid_fs) > 1
    text(0.1, 0.6, sprintf('FS统计 (帧):'), 'FontSize', 11, 'FontWeight', 'bold');
    text(0.1, 0.5, sprintf('  平均: %.1f', mean(valid_fs)), 'FontSize', 10);
    text(0.1, 0.4, sprintf('  标准差: %.1f', std(valid_fs)), 'FontSize', 10);
end

if length(valid_fo) > 1
    text(0.1, 0.3, sprintf('FO统计 (帧):'), 'FontSize', 11, 'FontWeight', 'bold');
    text(0.1, 0.2, sprintf('  平均: %.1f', mean(valid_fo)), 'FontSize', 10);
    text(0.1, 0.1, sprintf('  标准差: %.1f', std(valid_fo)), 'FontSize', 10);
end

xlim([0, 1]); ylim([0, 1]);
title('检测结果统计');
axis off;

fprintf('\n=== 分析完成 ===\n');
if length(valid_fs) > 0 || length(valid_fo) > 0
    fprintf('成功检测到步态事件！检查图形窗口查看详细结果。\n');
else
    fprintf('所有算法都未能检测到步态事件，可能需要调整参数或检查数据质量。\n');
end

%% 方法1: 基本调用（最简单）
fprintf('\n=== 创建基本可视化图 ===\n');

% 准备原始数据结构体（添加到您现有的代码中）
raw_markers = struct();
raw_markers.(heelMarkerName) = Markers_Corrected.(heelMarkerName);  % 使用滤波前数据
raw_markers.(toeMarkerName) = Markers_Corrected.(toeMarkerName);
raw_markers.SACR = Markers_Corrected.SACR;

% 准备滤波后数据结构体
filtered_markers = struct();
filtered_markers.(heelMarkerName) = filtheelmarker;
filtered_markers.(toeMarkerName) = filttoemarker;
filtered_markers.SACR = filtsacrmarker;

% 基本调用
fig_handles = plotGaitAnalysis('RawData', raw_markers, ...
                              'FilteredData', filtered_markers, ...
                              'SamplingFreq', f, ...
                              'Side', side, ...
                              'WalkingSpeed', vel2, ...
                              'HeelMarkerName', heelMarkerName, ...
                              'ToeMarkerName', toeMarkerName);