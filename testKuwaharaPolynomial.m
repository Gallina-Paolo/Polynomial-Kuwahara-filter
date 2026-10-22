function testKuwaharaPolynomial
%TESTKUWAHARAPOLYNOMIAL Check derivatives, planar equivalence and rank failures.
    [x,y] = meshgrid(0:14,0:12);
    I = 0.2+0.002*x+0.003*y+0.0004*x.^2+0.0002*x.*y+0.0003*y.^2;
    [J,gx,gy,info] = kuwaharaPolynomial(I,3,2);
    rr = 4:10; cc = 4:12;
    assertClose(J(rr,cc),I(rr,cc));
    assertClose(gx(rr,cc),0.002+0.0008*x(rr,cc)+0.0002*y(rr,cc));
    assertClose(gy(rr,cc),0.003+0.0002*x(rr,cc)+0.0006*y(rr,cc));
    assertClose(info.Derivatives(rr,cc,4),0.0008);
    assertClose(info.Derivatives(rr,cc,5),0.0002);
    assertClose(info.Derivatives(rr,cc,6),0.0006);
    rng(7); I = rand(9,11);
    [a,b,c,d] = kuwaharaPlanar(I,3);
    [e,f,g,h] = kuwaharaPolynomial(I,3,1);
    assertClose(a,e); assertClose(b,f); assertClose(c,g);
    assertClose(d.MinimumResidualVariance,h.MinimumResidualVariance);
    assert(isequal(d.Sector,h.Sector));
    [J,gx,gy] = kuwaharaPolynomial(0.4*ones(5),2,0);
    assertClose(J,0.4); assertClose(gx,0); assertClose(gy,0);
    assertFails(@() kuwaharaPolynomial(I,3,4),'kuwaharaPolynomial:Rank');
    assertFails(@() kuwaharaPolynomial(I,1,2),'kuwaharaPolynomial:Samples');
    I = 0.2+0.0001*x.^2.*y;
    [J,~,~,info] = kuwaharaPolynomial(I,3,3);
    assertClose(J(rr,cc),I(rr,cc));
    t = find(all(info.Exponents == [2 1],2));
    assertClose(info.Derivatives(rr,cc,t),0.0002);
    fprintf('All polynomial Kuwahara checks passed.\n');
end

function assertClose(a,b)
    delta = abs(a-b);
    assert(all(delta(:)<1e-10),'Unexpected numerical result.');
end

function assertFails(f,id)
    try
        f();
    catch err
        assert(strcmp(err.identifier,id),'Unexpected error: %s',err.identifier);
        return
    end
    error('An invalid input was accepted.');
end
