% Fast Poisson Editing
% Reference: http://www.cs.jhu.edu/~misha/Fall07/Papers/Perez03.pdf
% Based on Ramesh Raskar's script, available here:
% http://web.media.mit.edu/~raskar/photo/code.pdf

function poisson_reconstruct()
    % example usage:
    im_src = imread('i1.png');
    im_dst = imread('i2.png');
    sz = [700,700];
    im_src = imresize(im_src,sz);
    im_dst = imresize(im_dst,sz);
    im_res = blit_images(im_src,im_dst);

    clf();
    subplot(1,3,1);
    imshow(im_src);
    subplot(1,3,2);
    imshow(im_dst);
    subplot(1,3,3);  
    imshow(im_res);
end

function [Dx,Dy]=get_grads(im)
    % return the x and y gradients.
    imsz = size(im);
    H = imsz(1);
    W = imsz(2);
    Dx = zeros(H,W); 
    Dy = zeros(H,W);
    j = 1:H-1;
    k = 1:W-1;
    Dx(j,k) = im(j,k+1) - im(j,k);
    Dy(j,k) = im(j+1,k) - im(j,k);
end

function L = get_laplacian(Dx,Dy)
    % return the laplacian
    dsz = size(Dx);
    H = dsz(1);
    W = dsz(2);
    Dxx = zeros(H,W);
    Dyy = zeros(H,W);
    j = 1:H-1;
    k = 1:W-1;
    Dxx(j,k+1) = Dx(j,k+1) - Dx(j,k);
    Dyy(j+1,k) = Dy(j+1,k) - Dy(j,k);
    L = Dxx + Dyy;
end

function img=poisson_solve(Dx,Dy,bnd)
    bndsz = size(bnd);
    H = bndsz(1);
    W = bndsz(2);
    L = get_laplacian(Dx,Dy);

    % boundary image contains image intensities at boundaries
    % set the interior of the boundary-image to 0:
    bnd(2:end-1,2:end-1) = 0;
    j = 2:H-1;
    k = 2:W-1;
    L_bp = zeros(H,W);
    L_bp(j,k) = -4*bnd(j,k) + bnd(j,k+1) + bnd(j,k-1) + bnd(j-1,k) + bnd(j+1,k);
    clear j k
    L = L - L_bp;
    L = L(2:end-1,2:end-1);
    % compute the 2D DST:
    L_dst = dst(dst(L)')'; % first along the columns, then along rows
    % normalize:
    % compute Eigen Values
    [xx,yy] = meshgrid(1:W-2,1:H-2);
    D = (2*cos(pi*xx/(W-1))-2) + (2*cos(pi*yy/(H-1))-2);
    L_dst = L_dst ./ D;

    % put solution in inner points; outer points obtained from boundary image:
    img_interior = idst(idst(L_dst)')';
    img = bnd;
    img(2:end-1,2:end-1) = img_interior;
end

function x = clip(x,l,h) 
    x(x<l) = l;
    x(x>h) = h;
end

function im_res=blit_images(im_top,im_back,grad_shrink)
    %combine images using poission editing.
    %IM_TOP and IM_BACK should be of the same size.
    assert(all(size(im_top)==size(im_back)), 'the TOP and BACK images should be of the same size');
    im_top = double(im_top);
    im_back = double(im_back);

    im_res = zeros(size(im_top),'double');
    % process each channel independently:
    for ich = 1:size(im_top,3)
        ims = im_top(:,:,ich);
        imd = im_back(:,:,ich);
        % get the gradients of the two images:
        [gxs,gys] = get_grads(ims);
        [gxd,gyd] = get_grads(imd);

        % mix the source and target gradients:
        % this is based on the Perez et al.:
        gx = gxs;
        gxm = abs(gxd)>abs(gxs);
        gx(gxm) = gxd(gxm);
        gy = gys;
        gym = abs(gyd)>abs(gys);
        gy(gym) = gyd(gym);
        im_res(:,:,ich) = clip(poisson_solve(gx,gy,imd),0,255);
    end
    im_res = uint8(im_res);
end
