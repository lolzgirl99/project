newPackage( "FastMatrixEval",
Version => "0.1", Date => "February 20, 2019", Authors => {
     {Name => "Karl Schwede",
     Email=> "kschwede@gmail.com",
     HomePage=> "http://www.math.utah.edu/~schwede"
     },
     {Name => "Yuhui (Wei) Yao",
     Email=> "yuhuiyao4ever@gmail.com",
     }--
}, --this file is in the public domain
Headline => "A package for bertini and fast matrix evaluation.", DebuggingMode => true, Reload=>true)
export{
	"genericHyperplanePoint",
	"genericRankByHyperplane",
    "genericHyperplaneCurve",
    "embedCurveIn3Dim",
    "fastMinorsTableP",
    "fastMinorsTable",
	"Dim",
    "myRandom",
    "integralScheme",
    "Threads",
    "collapseMatrixV1",
    "collapseMatrixV2",
    --"TT1","TT2","TT3","TT4"
    "singDimMinor",
    "chooseSubmatrixLargestDegree",
    "MaxMinors",--an option, should be renamed
    "chooseMinorLargestDegree"
    "replaceZeros"
};
exportMutable{
    "TT1","TT2","TT3","TT4"
}

----------------------------------------------------------------
--************************************************************--
-------------------- Improving random --------------------------
--************************************************************--
----------------------------------------------------------------
fastMinors = (n, M1) -> (
    if (numColumns M1 > numRows M1) then M1 = transpose M1;
    ideal values fastMinorsTableP(n, M1)
);



temp = (Lis,n, M1, R1) -> (
    --print "temp started";
	K := apply(Lis, i -> (i => sum(0..n-1, curCol -> (-1)^(curCol)*M1_(i#0#(n-1), i#1#curCol)*R1#(toList drop(i#0, {n-1,n-1}), toList drop(i#1, {curCol, curCol})))));
	return K
);

fastMinorsTableP = method(Options=>{Threads=>(allowableThreads-2)});

fastMinorsTableP(ZZ, Matrix) := opts -> (n, M1) -> (
    if (allowableThreads <= 0) then return fastMinorsTable(n, M1);
    rowCt := numRows M1;
    colCt := numColumns M1;
    gg := subsets(rowCt, n);
    ff := subsets(colCt, n);
    L := toList ((set gg)**(set ff));
    lenL := length L;
    R1 := {};
    if (n > 2) then (R1 = fastMinorsTableP(n-1, submatrix'(M1, {rowCt-1}, ),opts)) else ( return new HashTable from  apply(L, i -> (i => det (M1^(i#0)_(i#1))) ) );
    tempL := null;
    --RList := {R1, new HashTable from R1};
    --L1 := take(L, {0,(length L-1)//4});
    --T1 := createTask(temp, (L1,n, M1, R1));
    --L2 := take(L, {(length L-1)//4, (length L-1)//2});
    --T2 := createTask(temp, (L2,n, M1, R1));
    --L3 := take(L, {(length L-1)//2, (3*(length L-1))//4});
    --T3 := createTask(temp, (L3,n, M1, R1));
    --L4 := take(L, {(3*(length L-1))//4, (length L)-1});
    --T4 := createTask(temp, (L4,n, M1, R1));
    --taskList := {};
    taskList := apply(opts.Threads, i -> (tempL =  take(L, {(i*(lenL-1))//(opts.Threads), ((i+1)*(lenL-1))//(opts.Threads)});
        return createTask(temp, (tempL, n, M1, R1) ); ) );
    --print "we are here";
    --schedule T1;
    --schedule T2;
    --schedule T3;
    --schedule T4;
    apply(taskList, t -> schedule t);
    while true do (--sleep 1; --print "wait";--print "i am waiting";
        --if isReady(T1) or isReady(T2) then print "something";
        if (all(taskList, t->isReady(t))) then break;
        --sleep 10;
        );
    --wait T1;
    --wait T2;
    --print "hello";
    --myList := flatten{taskResult(T1), taskResult(T2)};
    myList := flatten flatten apply(taskList, t -> taskResult(t));
    H := new HashTable from myList;
    --print "goodbye";
    --collectGarbage();
    apply(taskList, tt -> cancelTask(tt));
    return H;
)

fastMinorsTable = (n, M1) -> (
    rowCt := numRows M1;
    colCt := numColumns M1;
    gg := subsets(rowCt, n);
    ff := subsets(colCt, n);
    L := toList ((set gg)**(set ff));
    R := {};
    if (n > 2) then (R = fastMinorsTable(n-1, submatrix'(M1, {rowCt-1}, ))) else ( return new HashTable from  apply(L, i -> (i => det (M1^(i#0)_(i#1))) ) );
    H := new HashTable from  apply(L, i -> (i => sum(0..n-1, curCol -> (-1)^(curCol)*M1_(i#0#(n-1), i#1#curCol)*R#(toList drop(i#0, {n-1,n-1}), toList drop(i#1, {curCol, curCol}))))); --det (M^(i#0)_(i#1))

    return H;
)

myRandom = method(Options=>(options(random)++{Homogeneous=>true}));

myRandom(List, Ring) := RingElement => opts -> (L1, KR) -> (
    if (#L1 == 1) then (
        return myRandom(L1#0, KR, opts);
    ) else if (opts.Homogeneous == false) then ( --currently, we don't support random multidegrees
        error "myRandom: Expected a degree of length one for non-Homogeneous elements.";
    ) else (
        --the options for random do not include Homogeneous.  This mutable hashTable drops the Homogeneous option
        newOpts := new MutableHashTable from opts;
        remove(newOpts, Homogeneous);
        return random(L1, KR, new OptionTable from newOpts);
    )
);

myRandom(ZZ, Ring) := RingElement => opts -> (n1, KR) -> (
    --the options for random do not include Homogeneous.  This mutable hashTable drops the Homogeneous option
    newOpts := new MutableHashTable from opts;
    remove(newOpts, Homogeneous);
    if (opts.Homogeneous == true) then (
        return random(n1, KR, new OptionTable from newOpts);
        )
    else if ( n1 >= 0) then (
        finOpts := new OptionTable from newOpts;
        return sum(n1+1, i -> random(i, KR, finOpts));
        )
    else(
        error "myRandom: Expected a positive degree for non-Homogeneous elements.";
        );
);


myRandom(ZZ, FractionField) := RingElement => opts -> (n1, KR) -> (
   ambientRing := ring numerator sub(1, KR);
   num := myRandom(n1, ambientRing, opts);
   denom := myRandom(n1, ambientRing, opts);
   return sub(num/denom, KR);
);

myRandom(FractionField) := RingElement => opts -> (KR) -> (
   ambientRing := ring numerator sub(1, KR);
   num := random(opts.Height, ambientRing, opts);
   denom := random(opts.Height, ambientRing, opts);
   return sub(num/denom, KR);
);



--random := myRandom;


----------------------------------------------------------------
--************************************************************--
-------------------- Function Defitions ------------------------
--************************************************************--
----------------------------------------------------------------

genericHyperplanePoint = method(Options=>{Homogeneous=>false});

genericHyperplanePoint(Ring) := opts->R1 -> ( --
    ourDim := dim R1;
    genList := gens ambient R1;
    KK := coefficientRing(R1);
    ourIdeal := ideal R1;
    if (char KK > 0) then (
        --todo, figure out a good way to make kk bigger (how much bigger?)
    );
    S := ambient(R1);
        --todo, choose a general/random degree 1 element F of KK[genList] = S.
        O := random ourDim;
        RandomVar := genList#O;
        -- this picks a general variable without coefficient, need random element of KK
    if (opts.Homogeneous==true) then (
        --fi := null;
      --count := 0; i:=0; while (i < ourDim-1) do (i=i+1; fi = ideal(myRandom(1, S)); ourIdeal = ourIdeal+fi);
      ourIdeal = ourIdeal + sum(ourDim-1, i->ideal(myRandom(1, S, Homogeneous=>true)));
    )
    else (
      ourIdeal = ourIdeal + sum(ourDim, i->ideal(myRandom(1, S, Homogeneous=>false)));
    );
    return (S/ourIdeal);
);

--maybe we should create a function that creates the linear term, that the above and below just themselves use.
linearTerm = method(Options=>{});

linearTerm(Ring) := opts->R1 -> (
  ourDim := numgens ambient R1;
  genList := gens ambient R1;
  kk := coefficientRing(R1);
  ourIdeal := ideal R1;
  S := kk[genList];
  listElem := terms random(1, S);
  O := random ourDim;
  return listElem#O
);

genericRankByHyperplane = method(Options=>{Homogeneous=>false, Dim=>null, Prune=>true});

genericRankByHyperplane(Ring, Matrix) := opts -> (R1, M1) -> ( --
    --R1 := ring M1;
    d1 := opts.Dim;

    if (d1 === null) then (
        d1 = dim R1;
    );
    RH := trim genericHyperplanePoint(R1, Homogeneous=>opts.Homogeneous);
    phi := null;
    if (opts.Prune == true) then (prune RH; phi = (RH.minimalPresentationMap)*(map(RH, R1))) else (phi = map(RH, R1)); --if prune
    --print phi(M1);
    --1/0;
    return rank(phi(M1));
);

integralScheme = method(Options=>{Homogeneous=>false, Dim=>null, Prune=>true});

integralScheme(Ring) := opts -> (R1) -> (
  M1 := jacobian(ideal(R1));
  genRk := genericRankByHyperplane(R1, sub(M1, R1));
    if (genRk== numgens(ambient(R1))-dim(R1)  ) then (
      return true)
    else if (genRk<numgens(ambient R1)-dim(R1)  ) then return false;
);

genericHyperplaneCurve = method(Options=>{Homogeneous=>false});
genericHyperplaneCurve(Ring) := opts->R1 -> ( --
    ourDim := dim R1;
    genList := gens ambient R1;
    KK := coefficientRing(R1);
    ourIdeal := ideal R1;
    if (char KK > 0) then (
        --todo, figure out a good way to make kk bigger (how much bigger?)
    );
    S := ambient(R1);
        --todo, choose a general/random degree 1 element F of KK[genList] = S.
        O := random ourDim;
        RandomVar := genList#O;
        -- this picks a general variable without coefficient, need random element of KK
    if (opts.Homogeneous==true) then (
        --fi := null;
      --count := 0; i:=0; while (i < ourDim-1) do (i=i+1; fi = ideal(myRandom(1, S)); ourIdeal = ourIdeal+fi);
      ourIdeal = ourIdeal + sum(ourDim-2, i->ideal(myRandom(1, S, Homogeneous=>true)));
    )
    else (
      ourIdeal = ourIdeal + sum(ourDim-1, i->ideal(myRandom(1, S, Homogeneous=>false)));
    );
    return (S/ourIdeal);
);

embedCurveIn3Dim = method(Options=>{Homogeneous=>false});

embedCurveIn3Dim(Ring) := opts -> (R1) -> ( --you should pass this a 1-dimensional ring (or 2 dimensional, if homogeneous)
    coeffRing := coefficientRing R1;
   -- myMon := new monoid from [TT1]
    if (opts.Homogeneous == false) then (
        --myMon := new monoid from {TT1,TT2,TT3};
        S1 := coeffRing[TT1,TT2,TT3];
        phi1 := map(R1, S1, {myRandom(1, R1, Homogeneous=>false), myRandom(1, R1, Homogeneous=>false), myRandom(1, R1, Homogeneous=>false)});
        J1 := ker phi1; --this is the slow part
        return (dim (J1 + minors(2, jacobian J1)) < 0);
    )
    else if (opts.Homogeneous == true) then(
        S2 := coeffRing[TT1, TT2, TT3, TT4];
        phi2 := map(R1, S2, {myRandom(1, R1, Homogeneous=>true), myRandom(1, R1, Homogeneous=>true), myRandom(1, R1, Homogeneous=>true), myRandom(1, R1, Homogeneous=>true)});
        J2 := ker phi2;  --this is the slow part
        return (dim (J2 + minors(2, jacobian J2)) < 1);
    );
);

collapseMatrixV1 = method(Options=>{});

collapseMatrixV1(ZZ, Matrix) := opts -> (n1, M1) -> (--this randomly projects a matrix down to a smaller matrix
    R1 := ring M1;
    coeffRing := coefficientRing R1;
    transposeCount := 0;
    colList := apply(numColumns M1, i -> M1_i);
    --now lets make random columns
    randomCols := apply(n1, j->sum(#colList, i->(matrix (random(coeffRing)*(colList#i))) ));
    randomMat1 := fold((m1,m2)->m1|m2, randomCols);
    --print randomMat1;
    rowList := entries randomMat1;
    randomRows := apply(n1, j->sum(#rowList, i->(random(coeffRing)*matrix{rowList#i})));
    randomMat2 := fold((m1,m2)->m1||m2, randomRows);
    return randomMat2;
    --rowList := apply(numRows M1, i-> )
);

collapseMatrixV2 = method(Options=>{});

collapseMatrixV2(ZZ, Matrix) := opts -> (n1, M1) -> (--this randomly projects a matrix down to a smaller matrix
    R1 := ring M1;
    coeffRing := coefficientRing R1;
    transposeCount := 0;
    columnList := apply(numColumns M1, i -> M1_i);
    newRowList := null;
    newRowVec := null;
    newColumnList := null;
    newVec := null;
    a1 := numRows M1;
    b1 := numColumns M1;
    a2 := a1 - n1;
    b2 := b1 - n1;
    -- these are how many rows and columns respectively need to be cut down. The idea is now
    -- to do a pairwise sum.

    -- So, I'm not sure how to make k not equal j...
    count := 0;
    i:=0;
    k:=0;
    j:=0;
    while (i < b2) do (
        k = random (b1-i);
        newColumnList = drop(columnList, {k, k});
        j = random (b1-1-i);
        newVec =(random(coeffRing)*(columnList#k))+(newColumnList#j);
        columnList = insert(0, newVec, drop(newColumnList, {j, j}));
        i=i+1;
    );
    -- furthermore, is there a way to make "each" j, k different from one another in every loop?

    -- the following should be the concatonated matrix wrt columns:
    randomMatrix1 := fold((m1,m2)->(matrix m1)|(matrix m2), columnList);
    rowList2 := entries randomMatrix1;
    -- again how to make j=/= k
    count = 0; i=0;
    while (i < a2) do (
        k = random (a1-i);
        newRowList = drop(rowList2, {k, k});
        j = random (a1-1-i);
        newRowVec =(random(coeffRing)*(rowList2#k))+(newRowList#j);
        rowList2 = insert(0, newRowVec, drop(newRowList, {j, j}));
        i=i+1;
    );
    randomMatrix2 := matrix rowList2;--:= fold((m1,m2)->(matrix {m1})||(matrix {m2}), rowList2);
    return randomMatrix2;
);


chooseSubmatrixLargestDegree = method(Options => {});

chooseSubmatrixLargestDegree(ZZ, Matrix) := opts -> (n1, M1) -> (
    rCt := numRows M1;
    cCt := numColumns M1;
    i := 0;
    j := 0;
    curM1 := M1;
    --degreeMatrx := matrix apply(entries M1, l1 -> apply(l1, i->degree i));
    curList := null;
    curMax := null;
    curCol := null;
    curRow := null;
    newCurCol := null;
    newCurRow := null;

    returnRowList := {};--maybe these should be mutable lists, we'll do that later...
    returnColList := {};
    --keepRowList := {};
    --keepColList := {};

    while (i < n1) do (
        print concatenate("in loop, i =", toString(i));
        curList = flatten entries curM1;
        curMax = maxPosition curList;
        --print (curList#curMax);
        curRow = curMax // cCt;
        curCol = curMax % cCt;
        curM1 = submatrix'(curM1, {curRow}, {curCol});
        rCt = rCt - 1;
        cCt = cCt - 1;
        --we need to adjust based on submatrices
        newCurRow = curRow;
        --print curRow;
        j = (#returnRowList) - 1;
        --while (j >= 0) do (
            --if (newCurRow >= returnRowList#j) then (newCurRow = newCurRow + 1; print concatenate("j:", toString(j), " i:", toString(i)));
            --j = j-1;
        --);
        --and for columns
        newCurCol = curCol;
        j = (#returnColList) - 1;
        --while (j >= 0) do (
          --  if (newCurCol >= returnColList#j) then newCurCol = newCurCol + 1;
            --j = j-1;
        --);
        --keepRowList = append(returnRowList, curRow);
        returnRowList = append(returnRowList, newCurRow);
        returnColList = append(returnColList, newCurCol);
        i = i + 1;
    );
    --print returnRowList;
    returnRowList = new MutableList from returnRowList;
    returnColList = new MutableList from returnColList;
    j = n1-1;
    --print concatenate("j=",toString(j));
    while (j >= 0) do ( --compare entry j
        newCurRow = returnRowList#j;
        i = j-1;
        --print "here1";
        while  (i >= 0) do (--compare entry j to previous terms in the list
            --print "here2";
            if (newCurRow >= returnRowList#i) then (newCurRow = newCurRow + 1;);
            i = i-1;
        );
        returnRowList#j = newCurRow;
        j = j-1;
    );
    j = n1-1;
    while (j >= 0) do ( --compare entry j
        newCurCol = returnColList#j;
        i = j-1;
        while  (i >= 0) do (--compare entry j to previous terms in the list
            if (newCurCol >= returnColList#i) then (newCurCol = newCurCol + 1;);
            i = i-1;
        );
        returnColList#j = newCurCol;
        j = j-1;
    );
    return {new List from returnRowList, new List from returnColList};
);


chooseMinorLargestDegree = method(Options => {});

chooseMinorLargestDegree(ZZ, Matrix) := opts -> (n1, M1) -> (
  bestDegree := chooseSubmatrixLargestDegree(n1, M1);
  minorRowList := bestDegree#0;
  minorColList := bestDegree#1;

  return (M1^minorRowList)_minorColList;
);

replaceZeros= method(Options=>{});

replaceZeros(Matrix):= Matrix => o->(M2) -> (
  Mute := mutableMatrix M2;
  m := numRows M2;
  n := numColumns M2;
  largeDeg := (degree max flatten entries M2)#0;
  largeGen := (last gens ambient ring M2)^largeDeg;

  i := 0;

  while (i<n*m) do (
    Row := i//n;
    Col := i % n;
    if ((flatten entries(M2_{Col}))#Row==0)
    then Mute_(Row, Col)=largeGen;
    i=i+1;
    );
unMute := matrix Mute;
return unMute;
);

chooseSubmatrixSmallestDegree = method(Options=>{});

chooseSubmatrixSmallestDegree(ZZ, Matrix) := o -> (n1, M3) -> (
      M1 := replaceZeros(M3);
      rCt := numRows M1;
      cCt := numColumns M1;
      i := 0;
      j := 0;
      curM1 := M1;
      --degreeMatrx := matrix apply(entries M1, l1 -> apply(l1, i->degree i));
      curList := null;
      curMax := null;
      curCol := null;
      curRow := null;
      newCurCol := null;
      newCurRow := null;

      returnRowList := {};--maybe these should be mutable lists, we'll do that later...
      returnColList := {};

      while (i < n1) do (
          curList = flatten entries curM1;
          curMax = minPosition curList;
          --print (curList#curMax);
          curRow = curMax // cCt;
          curCol = curMax % cCt;
          curM1 = submatrix'(curM1, {curRow}, {curCol});
          rCt = rCt - 1;
          cCt = cCt - 1;
          --we need to adjust based on submatrices
          newCurRow = curRow;
          j = (#returnRowList) - 1;
          while (j >= 0) do (
              if (newCurRow >= returnRowList#j) then newCurRow = newCurRow + 1;
              j = j-1;
          );
          --and for columns
          newCurCol = curCol;
          j = (#returnColList) - 1;
          while (j >= 0) do (
              if (newCurCol >= returnColList#j) then newCurCol = newCurCol + 1;
              j = j-1;
          );
          returnRowList = append(returnRowList, newCurRow);
          returnColList = append(returnColList, newCurCol);
          i = i + 1;
          );
      return {returnRowList, returnColList};
    );


chooseMinorSmallestDegree = method(Options => {});

chooseMinorSmallestDegree(ZZ, Matrix) := o -> (n1, M1) -> (
    bestDegree := chooseSubmatrixSmallestDegree(n1, M1);
    minorRowList := bestDegree#0;
    minorColList := bestDegree#1;

    return (M1^minorRowList)_minorColList;
  );


















singDimMinor= method(Options=>{MaxMinors => 10});

singDimMinor(Ring, Matrix) := opts -> (R1, M1) -> (
  -- for now want dim of singular locus to be at least 2?
    ambR := ambient R1;
    Id := ideal R1;
    n := numgens ambR;
    r := dim R1;
    fullRank := n-r;
    workingM1 := new MutableMatrix from M1;

    largestSubmatrixLocation := chooseSubmatrixLargestDegree(fullRank, M1);
    rows1 := largestSubmatrixLocation#0;
    columns1 := largestSubmatrixLocation#1;

    minor1 := ideal values fastMinorsTable(fullRank, (M1^rows1)_columns1);
      -- take determinant of minor made using chooseLargestDegree
    quotient1 := ambR/(Id+sub(minor1, ambR));
    d := dim(quotient1);

    --largestSubmatrixLocation := chooseSubmatrixLargestDegree(fullRank, M1);
    biggestRow := rows1#0;
    biggestCol := columns1#0;
    --mNew := submatrix'(M1, rows1, columns1);
    workingM1_(biggestRow, biggestCol) = 0;
      -- dropping the rows columns used to make minor1, make new matrix

    i := 0;
    rowsNew := null;
    columnsNew := null;

    while ((r-d<2) and (i < opts.MaxMinors)) do (
      largestSubmatrixLocation = chooseSubmatrixLargestDegree(fullRank, new Matrix from workingM1);
      rows1 = largestSubmatrixLocation#0;
      columns1 = largestSubmatrixLocation#1;
      biggestRow = rows1#0;
      biggestCol = columns1#0;
      --mNew := submatrix'(M1, rows1, columns1);
      workingM1_(biggestRow, biggestCol) = 0;
      minor1 = minor1 + (ideal values fastMinorsTable(fullRank, (M1^rows1)_columns1) );
      -- make a new minor ideal out of the matrix after the rows/columns from the previous one are dropped
      -- add it to the old one
      -- take determinant
      quotient1 = ambR/(Id+sub(minor1, ambR));
      d = dim(quotient1);

    --rowsNew = chooseSubmatrixLargestDegree(fullRank, mNew)#1;
    --columnsNew = chooseSubmatrixLargestDegree(fullRank, mNew)#2;
    --mNew = submatrix'(mNew, rowsNew, columnsNew);
    -- make another submatrix after dropping rows/cols
      i = i+1;
    );
    1/0;
  --if (r-d>1) then (return "singular locus small");
  --if (r-d<2) then (return "singular locus large");
  --return {r, d};
);

--Can we do some Noether normalization?

--functionality to add

--1.
--figure out the (probabilistic) degree of a rational map
--choose several points randomly (by hyperplane sections) on Y.
--then given phi : X -> Y, let's find degree of inverse images of those points....

--2.
--Figure out stuff about the dimension of singular locus
--in particular, we should be able to detect whether a variety is R1, and thus we should be able to detect normality, at least probabilistically

--3.
--compare this with the speed of the point (rationalPoint) functionality in Cremona.  (And also compare accuracy).

--4.
--Implement our own rational point Noether normalization over bigger fields (find a polynomial subring such that the extension is a finite)
