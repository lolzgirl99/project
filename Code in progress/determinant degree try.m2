
determinantDeg=method(Options=>{});
determinantDeg(Matrix) := o->(M1) ->  (
  n1 := min(numColumns M1, numRows M1);
  l := apply(sort flatten entries M1, degree);
  l = select(l, i -> i>-infinity);
  return first apply(pack(n1, l), sum);
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

    while (i < n1) do (
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


newSingDimMinor= method(Options=>{MaxMinors => 100});

newSingDimMinor(Ring, Matrix) := opts -> (R1, M1) -> (
        ambR := ambient R1;
        Id := ideal R1;
        n := numgens ambR;
        r := dim R1;
        fullRank := n-r;
        myRand := 0;
        minTerm := sub(0, R1);
        colsOmit := numColumns M1- fullRank;
        mutM1 := new MutableMatrix from M1;
        mutM2 := new MutableMatrix from M1;

        submatrixL1 := chooseSubmatrixLargestDegree(fullRank, M1);
        mRowListL := submatrixL1#0;
        mColListL := submatrixL1#1;
        largestSubmatrix :=(M1^mRowListL)_mColListL;
        termList := flatten entries monomials(mutM1_(mRowListL#0, mColListL#0));
        if (#termList > 1) then ( --drop the smallest degree term basically, and sum what's left
            minTerm = minPosition(termList);
            mutM1_(mRowListL#0, mColListL#0) = sum(drop(termList, {minTerm, minTerm}));
        )
        else if (#termList == 1) then (
            minTerm = new MutableList from factor(termList#0);
            myRand = random(#minTerm);
            --myExp = (minTerm#myRand#1) - 1;
            minTerm#myRand = new Power from {(minTerm#myRand)#0, (minTerm#myRand#1) - 1};
            mutM1_(mRowListL#0, mColListL#0) = value (new Product from minTerm);
        )
        else(
            mutM1_(mRowListL#0, mColListL#0) = 0;
        );

        submatrixS1 := chooseSubmatrixSmallestDegree(fullRank, M1);
        mRowListS := submatrixS1#0;
        mColListS := submatrixS1#1;
        smallestSubmatrix :=(M1^mRowListS)_mColListS;
        mutM2_(mRowListS#0, mColListS#0) = (mutM2_(mRowListS#0, mColListS#0))*((gens ambR)#(random(#(gens ambR))));
            --multiply our mutable matrix by a random generator, for the future!

        minorL1 := determinant largestSubmatrix;
        minorS1 := determinant smallestSubmatrix;
        sumMinors := ideal(sub(minorL1, ambR))+ideal(sub(minorS1, ambR));
          -- take determinant
        quotient1 := ambR/(Id+sumMinors);
        d := dim(quotient1);

        i := 0;
        j := 1;
        minorL := null;
        minorS := null;
        nextMatrixL := M1;
        nextMatrixS := M1;

        while ((r-d<2) and (i < opts.MaxMinors)) do (
          while (i <= 2^j) do (
            --nextMatrixL = submatrix'(nextMatrixL, {(submatrixL1#0)#i},{(submatrixL1#1)#i});
            --nextMatrixS = submatrix'(nextMatrixS, {(submatrixS1#0)#i},{(submatrixS1#1)#i});

            submatrixL1 = chooseSubmatrixLargestDegree(fullRank, (new Matrix from mutM1));
            submatrixS1 = chooseSubmatrixSmallestDegree(fullRank, (new Matrix from mutM2));

            mRowListL = submatrixL1#0;
            mColListL = submatrixL1#1;
            largestSubmatrix =(M1^mRowListL)_mColListL;
            mRowListS = submatrixS1#0;
            mColListS = submatrixS1#1;
            smallestSubmatrix =(M1^mRowListS)_mColListS;

            minorL = determinant largestSubmatrix;
            minorS = determinant smallestSubmatrix;
            sumMinors = sumMinors + sub(minorL, ambR) + sub(minorS, ambR);
            i=i+1;
            );
          quotient1 = ambR/(Id+sumMinors);
          d = dim(quotient1);
          j = j+1;
        );
        return d;
    );
