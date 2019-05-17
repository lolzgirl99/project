randomMatrix = method(Options=>{Homogeneous => true});

randomMatrix(ZZ, ZZ, ZZ, Ring) := (n, m, d, R) -> (

    return 7;
);


randomEval = method(Options=>{Height=>100, CoeffRing=>null});

randomEval(Matrix) := Matrix => opts->(M1) -> (
        R1 := ring M1;
        myCoeff := opts.CoeffRing;
        if (myCoeff === null) then (
            myCoeff = coefficientRing R1;
            );
        numVars := numColumns vars R1;
        randomMap := map(myCoeff, R, toList apply(numVars, i->random(myCoeff, Height=>opts.Height)));
        return randomMap(M1);
);

makeRandomMatrix = method(Options=>{Height=>10});

makeRandomMatrix(ZZ, ZZ, ZZ, Ring) := Matrix => opts -> (n, m, d, R1) -> (
    return sum(apply(d+1, j->random(R1^(apply(n, i->j)), R1^m, Height=>o.Height)));
    );

symbol Projection;
symbol Hyperplane;

findRandomPointOnAffine = method(Options=>{Strategy=>Projection, CoefficientRing=>null});

findRandomPointOnAffine(Ring) := opts -> (R1) -> (
    myCoeff := opts.CoefficientRing;
    if (myCoeff === null) then myCoeff = coefficientRing(R1);
    if (opts.Strategy == Projection) then (
        return findPtOnAffineByProjection(R1, CoefficientRing => myCoeff);
    )
    else if (opts.Strategy == Hyperplane) then (

    )
);

findPtOnAffineByHyperplane= method(Options=>{CoefficientRing=>null});

findPtOnAffineByHyperplane(Ring) := opts -> (R1) -> (
    
);

needsPackage "NoetherNormalization";
needsPackage "InvolutiveBases"; --for NoetherNormalization

findPtOnAffineByProjection = method(Options=>{CoefficientRing=>null});

findPtOnAffineByProjection(Ring) := opts -> (R1) -> (
    --NNList := noetherNormalization(R1);

    T1 := ambient R1;
    --NNList := invNoetherNormalization(ideal R1);
    --S1 := (opts.CoefficientRing)[NNList#1];
    --NNMap := (map(R1,T1))*(map(T1, T1, NNList#0))*(map(T1, S1));
    NNList := noetherNormalization(R1);
    S1 := (opts.CoefficientRing)[NNList#2];
    NNMap := (map(R1,T1))*(NNList#0)*(map(T1, S1)); --this is the map from the Noether
                                                                     --Normalization to to polynomial ring,
                                                                     --composed with the transform of the polynomial rings
                                                                     --composed witht he project to our actual ring
                                                                     --It should be the Noether Normalization map.
    --return NNMap;
    print NNMap;
    d := #(gens(S1));
    Ip := NNMap(ideal(apply(d, j -> random(1, S1) + random(0, S1))));
    return Ip;
);
