%% cross-trial alignment

A = double(Ding1);
B = double(Ding2);

nb_tiles = 16;
% [xtable, ytable, utable, vtable, typevector] = piv_DCC(A,B,8,256/nb_tiles,1,[],[]);
% [xtable, ytable, utable, vtable, typevector] = piv_FFTmulti(image1,image2,interrogationarea, step, subpixfinder, mask_inpt, roi_inpt,passes,int2,int3,int4,imdeform);
[xtable, ytable, utable, vtable, typevector] = piv_FFTmulti(A,B,32,256/nb_tiles, 1, [], [],2,32,64,128,'linear');

v_pixels{k} = smoothn(vtable,0.3);
u_pixels{k} = smoothn(utable,0.3);
mean(u_pixels{k}(:))
[X,Y] = meshgrid(1:size(A,2), 1:size(A,1));
v_interpol{k} = interp2( xtable,ytable,v_pixels{k},X,Y );
u_interpol{k} = interp2( xtable,ytable,u_pixels{k},X,Y );