
determinantDeg=method(Options=>{});
determinantDeg(Matrix) := o->(M1) ->  (
  n1 := min(numColumns M1, numRows M1);
  l := apply(sort flatten entries M1, degree);
  l = select(l, i -> i>-infinity);
  return first apply(pack(n1, l), sum);
);

chooseSubmatrixLargestDegree = method(Options => {});

chooseSubmatrixLargestDegree(ZZ, MutableMatrix) := opts -> (n1, M1) -> (
    rCt := numRows M1;
    cCt := numColumns M1;
    i := 0;
    j := 0;
    curM1 := new Matrix from M1;
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
        --print concatenate("in loop, i =", toString(i));
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
chooseMinorLargestDegree(ZZ, Matrix) := o -> (n1, M1) -> (
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
  largeDeg := sum((degree max flatten entries M2));
  --print largeDeg;
  largeGen := (last gens ambient ring M2)^(2*largeDeg+1);

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

chooseSubmatrixSmallestDegree(ZZ, MutableMatrix) := o -> (n1, M3) -> (
          --M1 := new Matrix replaceZeros(M3);
          M1 := new Matrix from M3;
          rCt := numRows M1;
          cCt := numColumns M1;
          i := 0;
          j := 0;
          curM1 := M1;
          curList := null;
          curMax := null;
          curCol := null;
          curRow := null;
          newCurCol := null;
          newCurRow := null;

          returnRowList := {};
          returnColList := {};

          while (i < n1) do (
              --print concatenate("in loop, i =", toString(i));
              curList = flatten entries curM1;
              curMax = minPosition curList;
              curRow = curMax // cCt;
              curCol = curMax % cCt;
              curM1 = submatrix'(curM1, {curRow}, {curCol});
              rCt = rCt - 1;
              cCt = cCt - 1;
              newCurRow = curRow;
              j = (#returnRowList) - 1;
              newCurCol = curCol;
              j = (#returnColList) - 1;
              returnRowList = append(returnRowList, newCurRow);
              returnColList = append(returnColList, newCurCol);
              i = i + 1;
          );
          returnRowList = new MutableList from returnRowList;
          returnColList = new MutableList from returnColList;
          j = n1-1;
          while (j >= 0) do (
              newCurRow = returnRowList#j;
              i = j-1;
              while  (i >= 0) do (
                  if (newCurRow >= returnRowList#i) then (newCurRow = newCurRow + 1;);
                  i = i-1;
              );
              returnRowList#j = newCurRow;
              j = j-1;
          );
          j = n1-1;
          while (j >= 0) do (
              newCurCol = returnColList#j;
              i = j-1;
              while  (i >= 0) do (
                  if (newCurCol >= returnColList#i) then (newCurCol = newCurCol + 1;);
                  i = i-1;
              );
              returnColList#j = newCurCol;
              j = j-1;
          );
          return {new List from returnRowList, new List from returnColList};
      );

chooseMinorSmallestDegree = method(Options => {});
chooseMinorSmallestDegree(ZZ, MutableMatrix) := o -> (n1, M1) -> (
    bestDegree := chooseSubmatrixSmallestDegree(n1, M1);
    minorRowList := bestDegree#0;
    minorColList := bestDegree#1;
    return (M1^minorRowList)_minorColList;
  );


replaceSmallestTerm= method(Options=>{});
replaceSmallestTerm(List, MutableMatrix) := opts -> (submatrixS, M1) -> (
  ambR:= ambient ring(M1);
  --submatrixS := chooseSubmatrixSmallestDegree(FR, M1);
  rowListS := submatrixS#0;
  colListS := submatrixS#1;
  mutedSM := new MutableMatrix from M1;
  moddedRow := random(#rowListS);
  moddedCol := random(#colListS);
  --mutedSM_(rowListS#0, colListS#0) = (mutedSM_(rowListS#0, colListS#0))*((gens ambR)#(random(#(gens ambR))));
  mutedSM_(moddedRow, moddedCol) = (mutedSM_(moddedRow, moddedCol))*((gens ambR)#(random(#(gens ambR))));
  --unMuteSM := new Matrix from mutedSM;
  return mutedSM;
  );

replaceLargestTerm= method(Options=>{});
replaceLargestTerm(List, MutableMatrix):= opts-> (submatrixL, M1)->(
  --submatrixL := chooseSubmatrixLargestDegree(FR, M1);
  rowListL := submatrixL#0;
  colListL := submatrixL#1;
  largeSubmatrix :=(M1^rowListL)_colListL;
  termLists := flatten entries monomials(M1_(rowListL#0, colListL#0));

  muted := new MutableMatrix from M1;
  if (#termLists > 1) then ( --drop the smallest degree term basically, and sum what's left
      minTerm = minPosition(termLists);
      muted_(rowListL#0, colListL#0) = sum(drop(termLists, {minTerm, minTerm}));
  )
  else if (#termLists == 1) then (
      minTerm = new MutableList from factor(termLists#0);
      if (#minTerm > 0) then (
          myRand = random(#minTerm);
          minTerm#myRand = new Power from {(minTerm#myRand)#0, (minTerm#myRand#1) - 1};
          muted_(rowListL#0, colListL#0) = value (new Product from minTerm);
      )
      else (
          muted_(rowListL#0, colListL#0) = sub(0, ring M1);
      )
  )
  else(
      muted_(rowListL#0, colListL#0) = 0;
      );
  --replaced := new Matrix from muted;
  return muted;
  );


newSingDimMinor= method(Options=>{MaxMinors => 1000});

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
        mutM2 := new MutableMatrix from replaceZeros(M1);

        submatrixL1 := chooseSubmatrixLargestDegree(fullRank, mutM1);
        mRowListL := submatrixL1#0;
        mColListL := submatrixL1#1;
        largestSubmatrix :=(M1^mRowListL)_mColListL;

        submatrixS1 := chooseSubmatrixSmallestDegree(fullRank, mutM2);
        mRowListS := submatrixS1#0;
        mColListS := submatrixS1#1;
        smallestSubmatrix :=(M1^mRowListS)_mColListS;

        minorL1 := determinant largestSubmatrix;
        minorS1 := determinant smallestSubmatrix;
        sumMinors := sub(minorL1, ambR)+sub(minorS1, ambR);
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
            mutM1 = replaceLargestTerm(submatrixL1, mutM1);
            mutM2 = replaceSmallestTerm(submatrixS1, mutM2);
            print submatrixS1;

            submatrixL1 = chooseSubmatrixLargestDegree(fullRank, mutM1);
            submatrixS1 = chooseSubmatrixSmallestDegree(fullRank, mutM2);

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
        print new Matrix from mutM1;
        print new Matrix from mutM2;
        return d;
    );
