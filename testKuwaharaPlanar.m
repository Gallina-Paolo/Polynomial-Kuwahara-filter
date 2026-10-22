function testKuwaharaPlanar
%TESTKUWAHARAPLANAR Check RGB conversion, affine reproduction and invalid inputs.
    rng(7);
    rgb = uint8(randi([0 255],12,15,3));
    d = double(rgb)/255;
    grayscale = 0.2989360213*d(:,:,1)+0.5870430745*d(:,:,2)+0.1140209042*d(:,:,3);
        [J,gx,gy] = kuwaharaPlanar(rgb,3);
        [J0,gx0,gy0] = kuwaharaPlanar(grayscale,3);
        assertClose(J,J0); assertClose(gx,gx0); assertClose(gy,gy0);
        assert(isequal(size(J),[12 15]));
    [x,y] = meshgrid(0:14,0:11);
    I = 0.2+0.01*x+0.02*y;
    [J,gx,gy] = kuwaharaPlanar(I,3);
    assertClose(J(4:end-3,4:end-3),I(4:end-3,4:end-3));
    assertClose(gx(4:end-3,4:end-3),0.01);
    assertClose(gy(4:end-3,4:end-3),0.02);
    assertClose(kuwaharaPlanar(repmat(I,1,1,3),3),kuwaharaPlanar(I,3));
    assertClose(kuwaharaPlanar(true(5),3),ones(5));
    assertClose(kuwaharaPlanar(uint16(65535*ones(5)),3),ones(5));
    for input = {[], zeros(3,3,4), NaN(3), 'foto.png'}
        failed = false;
        try
            kuwaharaPlanar(input{1},3);
        catch
            failed = true;
        end
        assert(failed,'An invalid input was accepted.');
    end
    fprintf('All planar Kuwahara checks passed.\n');
end

function assertClose(a,b)
    delta = abs(a-b);
    assert(all(delta(:)<1e-11),'Unexpected numerical result.');
end

