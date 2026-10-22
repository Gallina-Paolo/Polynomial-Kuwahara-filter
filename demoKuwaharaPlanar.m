%DEMOKUWAHARAPLANAR Standalone synthetic example; no external images required.
rng(7);
[x,y] = meshgrid(linspace(0,1,160), linspace(0,1,120));
ideal = 0.15 + 0.25*x + 0.15*y + 0.25*(x > 0.55);
I = min(1,max(0,ideal + 0.025*randn(size(ideal))));
radius = 3;
[J,gx,gy,info] = kuwaharaPlanar(I, radius);
figure('Name','Planar regression Kuwahara');
subplot(2,2,1); imagesc(I,[0 1]); axis image off; title('Input');
subplot(2,2,2); imagesc(J,[0 1]); axis image off; title('Filtered');
subplot(2,2,3); imagesc(info.Magnitude); axis image off; colorbar; title('Gradient magnitude');
subplot(2,2,4); imagesc(info.MinimumResidualVariance); axis image off; colorbar; title('Residual variance');
colormap(gray(256));
