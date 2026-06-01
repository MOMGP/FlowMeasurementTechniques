
clc, clear, close

window_size = 32;
pixel_to_mm = 1; %TODO - must change
delta_t = 0.01;
%TODO - change

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

        V(j,i,1) = dx *pixel_to_mm/delta_t;
        V(j,i,2) = dy * pixel_to_mm/delta_t;
    end
end

plot_velocity_field(V, windowed_hor_size, windowed_vert_size, window_size);

%TODOs still - add overlap, maybe multipass, FFT, wtf is normalized cross
%correlation??
%sub-pixel peak fitting???


%% Functions

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


function plot_velocity_field(V, windowed_hor_size, windowed_vert_size, window_size)

    [X,Y] = meshgrid( ...
        (1:windowed_hor_size)*window_size - window_size/2, ...
        (1:windowed_vert_size)*window_size - window_size/2);
    %disp(max(V(:))) %for debugging
    v_x = V(:,:,1);
    v_y = V(:,:,2);

    figure
    
    quiver(X, flipud(Y), v_x, flipud(v_y), 'r');

    axis equal
    drawnow
end


