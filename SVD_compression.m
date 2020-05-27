img = double(rgb2gray(imread('lena_Std.tif')));

[U,S,V]=svd(img);

figure;

for i = 1:size(img,1)

    imagesc(U(:,1:i)*S(1:i,1:i)*V(:,1:i)');

    colormap('gray')

    title(['layers added upto ', num2str(i)]);

    if i<30

        pause(0.2)

    elseif (30<=i)&&(i<100)

        pause(0.1)

    else

        pause(0.01)

    end

end