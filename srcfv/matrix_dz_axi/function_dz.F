
    distm1    = ONE/abs( yc(i,j) - ym )

    func0(i,j,1) = w(i,j,4)
    func0(i,j,2) = w(i,j,2) * velz(i,j) - mu(i,j) * gradw(i,j,1)
    func0(i,j,3) = w(i,j,3) * velz(i,j) - mu(i,j) * gradw(i,j,2) + mu(i,j) * velz(i,j) * distm1
    func0(i,j,4) = w(i,j,4) * velz(i,j) + p(i,j) + TWOTHIRD * mu(i,j) * ( gradu(i,j,1) + gradv(i,j,2) - TWO * vely(i,j) * distm1 )
    func0(i,j,5) = ( w(i,j,5) + p(i,j) ) * velz(i,j)

    func1(i,j,1) = ZERO
    func1(i,j,2) = velz(i,j)
    func1(i,j,3) = velz(i,j)
    func1(i,j,4) = velx(i,j)
    func1(i,j,5) = velz(i,j)

    func2(i,j,1) = ZERO
    func2(i,j,2) = gradw(i,j,1)
    func2(i,j,3) = gradw(i,j,2)
    func2(i,j,4) = vely(i,j)
    func2(i,j,5) = gradw(i,j,1)

    func3(i,j,1) = ZERO
    func3(i,j,2) = ZERO
    func3(i,j,3) = ZERO
    func3(i,j,4) = gradu(i,j,1) + gradv(i,j,2)
    func3(i,j,5) = velx(i,j)

    func4(i,j,1) = ZERO
    func4(i,j,2) = ZERO
    func4(i,j,3) = ZERO
    func4(i,j,4) = ZERO
    func4(i,j,5) = gradu(i,j,1)

    func5(i,j,1) = ZERO
    func5(i,j,2) = ZERO
    func5(i,j,3) = ZERO
    func5(i,j,4) = ZERO
    func5(i,j,5) = velz(i,j)

    func6(i,j,1) = ZERO
    func6(i,j,2) = ZERO
    func6(i,j,3) = ZERO
    func6(i,j,4) = ZERO
    func6(i,j,5) = gradw(i,j,2)

    func7(i,j,1) = ZERO
    func7(i,j,2) = ZERO
    func7(i,j,3) = ZERO
    func7(i,j,4) = ZERO
    func7(i,j,5) = vely(i,j)

    func8(i,j,1) = ZERO
    func8(i,j,2) = ZERO
    func8(i,j,3) = ZERO
    func8(i,j,4) = ZERO
    func8(i,j,5) = gradv(i,j,2)

    func9(i,j,1) = ZERO
    func9(i,j,2) = ZERO
    func9(i,j,3) = ZERO
    func9(i,j,4) = ZERO
    func9(i,j,5) = mu(i,j) * velx(i,j)

    func10(i,j,1) = ZERO
    func10(i,j,2) = ZERO
    func10(i,j,3) = ZERO
    func10(i,j,4) = ZERO
    func10(i,j,5) = mu(i,j) * vely(i,j)

    func11(i,j,1) = ZERO
    func11(i,j,2) = ZERO
    func11(i,j,3) = ZERO
    func11(i,j,4) = ZERO
    func11(i,j,5) = gradw(i,j,2) - velz(i,j) * distm1

    func12(i,j,1) = ZERO
    func12(i,j,2) = ZERO
    func12(i,j,3) = ZERO
    func12(i,j,4) = ZERO
    func12(i,j,5) = mu(i,j) * velz(i,j)

    func13(i,j,1) = ZERO
    func13(i,j,2) = ZERO
    func13(i,j,3) = ZERO
    func13(i,j,4) = ZERO
    func13(i,j,5) = gradu(i,j,1) + gradv(i,j,2) - TWO * vely(i,j) * distm1

    func14(i,j,1) = ZERO
    func14(i,j,2) = ZERO
    func14(i,j,3) = ZERO
    func14(i,j,4) = ZERO
    func14(i,j,5) = ZERO

    func15(i,j,1) = ZERO
    func15(i,j,2) = ZERO
    func15(i,j,3) = ZERO
    func15(i,j,4) = ZERO
    func15(i,j,5) = ZERO




