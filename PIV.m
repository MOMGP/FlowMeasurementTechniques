
clc, clear, close all

window_size = 32;
delta_t = 74 * 10^(-6);
%TODO - change

calibration_file = "data\PIV\FMT Results\cal_final\B00001.tif";
pixel_to_mm = calibration(calibration_file, false);
disp(["pix per mm", 1/pixel_to_mm])
[img1, img2] = split_image("data\PIV\FMT Results\AoA15_final\B00001.tif");
%figure(1)
%imshow(img1)
%figure(2)
%imshow(img2)
%figure(3)
%windowed_img = blockproc(img1, ...
%    [window_size window_size], ...
%    @(block) sum(block.data(:)));

%imshow(windowed_img,[])
[height, width] = size(img1);

windowed_vert_size = floor(height/window_size);
windowed_hor_size = floor(width/window_size);

V = zeros(windowed_vert_size, windowed_hor_size, 2);

for i = 1:windowed_hor_size
    for j = 1:windowed_vert_size
        %V(i,j) = find_displacement_of_window(img1, img2, window_size, j, i)/delta_t;
        [dx, dy] = find_displacement_of_window(img1, img2, window_size, i, j);
        %fprintf("dx=%f dy=%f\n", dx, dy); %for debugging

        V(j,i,1) = dx *pixel_to_mm/delta_t / 1000; %to convert to m/s
        V(j,i,2) = dy * pixel_to_mm/delta_t / 1000; %to convert to m/s
    end
end

plot_velocity_field(V, windowed_hor_size, windowed_vert_size, window_size, pixel_to_mm);

%TODOs still - add overlap, maybe multipass, FFT, wtf is normalized cross
%correlation??
%sub-pixel peak fitting???

plot_processed_data("data\PIV processed\Processed15\Overlap0SinglePass\B00001.dat");

%% Functions

function [pix_to_mm] = calibration(file_name, plot)
    img = imread(file_name);
    if plot
        figure
        imshow(img, []);
        hold on
    end
    
    x1 = 330;
    x2 = 1050;
    y1 = 1160;
    y2 = 1150;
    
    if plot
    plot([x1, x2], [y1, y2]) %set to 80 mm (I think)
    end
    pix_to_mm = 80/sqrt((x2-x1)^2 + (y2-y1)^2);
end


function [img1, img2] = split_image(file_name)
    img = imread(file_name);
    [height, width] = size(img);
    img1 = imcrop(img, [0, floor(height/2)+1, width, floor(height/2)]);
    img2 = imcrop(img, [0, 0, width, floor(height/2)]);
end
function [x_disp_pix, y_disp_pix] = find_displacement_of_window(img1, img2, window_size, wind_index_x, wind_index_y)
    x_pixels = (wind_index_x-1) * window_size + 1 : wind_index_x * window_size;
    y_pixels = (wind_index_y-1) * window_size +1 : wind_index_y * window_size;
    
    %for more numerical accuracy
    window1 = double(img1(y_pixels, x_pixels));
    window2 = double(img2(y_pixels, x_pixels));
    
    %maybe add normalization:
    window1 = double(window1) - mean(window1(:));
    window2 = double(window2) - mean(window2(:));
    
    corr_map = xcorr2(window1, window2); %todo - check if vectors are in reversed order. If so, change the order
    [~, vectorized_index] = max(corr_map(:));
    %ind_y = ceil(vectorized_index/(2*window_size-1));
    %ind_x = mod(vectorized_index, (2*window_size-1))+1;
    [ind_y, ind_x] = ind2sub(size(corr_map), vectorized_index);
    x_disp_pix = ind_x - window_size;
    y_disp_pix = ind_y - window_size;
end


function plot_velocity_field(V, windowed_hor_size, windowed_vert_size, window_size, pixel_to_mm, mask)
    if nargin < 6 || isempty(mask)
        mask = false(windowed_vert_size, windowed_hor_size);
    end
    [X,Y] = meshgrid( ...
        (1:windowed_hor_size)*window_size - window_size/2, ...
        (1:windowed_vert_size)*window_size - window_size/2);
    X = X * pixel_to_mm;
    Y = Y * pixel_to_mm;
    
    size(X)
    %disp(max(V(:))) %for debugging
    v_x = V(:,:,1);
    v_y = V(:,:,2);

    if mask
        % ensure mask is logical
        mask = logical(mask);
    
        % apply mask to velocity components
        v_x(~mask) = NaN;
        v_y(~mask) = NaN;
    end

    vel_magnitudes = flipud(hypot(v_x, v_y));
    
    
      

    figure
    %plotting the magnitude in terms of background color - like DaVis
    imagesc(X(1,:), Y(:,1), vel_magnitudes)
    set(gca,'YDir','normal')
    axis equal
    colormap(parula)
    colorbar
    caxis([0 20]);

    hold on
    quiver(X, flipud(Y), v_x, flipud(v_y), 'r');

    drawnow

end



function [] = plot_processed_data(file_name)
    %for all dat files, I hope NumHeaderLines=3
    data = importdata(file_name, ' ', 3);

    M = data.data;

    X   = M(:,1);
    Y   = M(:,2);
    v_x = M(:,3);
    v_y = M(:,4);
    
    
    % need to reshape cuz they are all flattened
    [n_cols, n_rows] = extract_cols_rows(file_name);

    X  = reshape(X, n_cols, n_rows)';
    Y  = reshape(Y, n_cols, n_rows)';
    v_x = reshape(v_x, n_cols, n_rows)';
    v_y = reshape(v_y, n_cols, n_rows)';



    % plotting like plot_vel_field
    vel_magnitudes = hypot(v_x, v_y);
    figure
    imagesc(X(1,:), Y(:,1), vel_magnitudes)
    set(gca,'YDir','normal')
    axis equal
    colormap(parula)
    colorbar

    hold on
    quiver(X, Y, v_x, v_y, 'r')

    drawnow
end

function [n_cols, n_rows] = extract_cols_rows(file_name)
    file_index = fopen(file_name,'r');
    
    fgetl(file_index); % title - first line
    fgetl(file_index); % vars - second line
    third_line = fgetl(file_index); % third line ( where it has I and J).

    fclose(file_index);
    % reading last header line using regex
    reg_expr = regexp(third_line, 'I=(\d+),\s*J=(\d+)', 'tokens');
    n_cols = str2double(reg_expr{1}{1}); % I
    n_rows = str2double(reg_expr{1}{2}); % J
end