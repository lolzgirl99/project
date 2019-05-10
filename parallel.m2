makeTableOf2x2Minors = M -> (
    rowCt = numRows M;
    colCt = numColumns M;
    i = 0;
    j = 0;
    k = 0;
    l = 0;

    L = flatten apply(rowCt-1, i -> flatten toList apply(i+1..rowCt-1, j-> flatten toList apply(colCt-1, k -> toList apply(k+1..colCt-1, l -> {{i,j},{k,l}}))) );
    H = new HashTable from apply(L, elt -> (elt => det (M^(elt#0)_(elt#1))));

    return H
)

makeTableOf3x3Minors = M -> (
    rowCt := numRows M;
    colCt := numColumns M;
  G = flatten toList apply(rowCt-1, i -> flatten toList apply (i+1..(rowCt-1), j -> flatten toList apply (j+1..(rowCt-1), k -> flatten  toList apply (colCt-1, l -> flatten toList apply(l+1..(colCt-1), m -> toList apply(m+1..(colCt-1), n -> {{i,j,k},{l,m,n}}))))));
    H = makeTableOf2x2Minors(M);
    T = flatten toList apply(rowCt-1, i -> flatten toList apply (i+1..(rowCt-1), j -> flatten toList apply (j+1..(rowCt-1), k -> flatten  toList apply (colCt-1, l -> flatten toList apply(l+1..colCt-1, m -> flatten toList apply(m+1..colCt-1, n -> {(i,l), {{j,k},{m,n}}, (i,m), {{j,k},{l,n}}, (i,n), {{j,k},{l,m}}} ))))));
    B = for i when i < ((length T)-5) list (if i%6 != 0 then continue; ((((M_(T#i))*H#(T#(i+1))) - ((M_(T#(i+2)))*H#(T#(i+3))) + ((M_(T#(i+4)))*H#(T#(i+5))))));
    P = new HashTable from apply((length G), i -> (G#i => B#i));

    return P
)

makeTableOf4x4Minors = M -> (
  H = makeTableOf3x3Minors(M);
  G = flatten toList apply(rowCt-1, i -> flatten toList apply (i+1..rowCt-1, j -> flatten toList apply (j+1..rowCt-1, k -> flatten toList apply (k+1..rowCt-1, l -> flatten  toList apply (colCt-1, m -> flatten toList apply(m+1..colCt-1, n -> flatten toList apply(n+1..colCt-1, o -> toList apply(o+1..colCt-1, p -> {{i,j,k,l},{m,n,o,p}}))))))));
  T = flatten toList apply(rowCt-1, i -> flatten toList apply (i+1..rowCt-1, j -> flatten toList apply (j+1..rowCt-1, k -> flatten toList apply (k+1..rowCt-1, l -> flatten  toList apply (colCt-1, m -> flatten toList apply(m+1..colCt-1, n -> flatten toList apply(n+1..colCt-1, o -> flatten toList apply(o+1..colCt-1, p -> {(i,m), {{j,k,l},{n,o,p}}, (i,n), {{j,k,l},{m,o,p}}, (i,o), {{j,k,l},{m,n,p}}, (i,p), {{j,k,l},{m,n,o}}}))))))));
  B = for i when i < ((length T)-7) list (if i%8 != 0 then continue; ((((M_(T#i))*H#(T#(i+1))) - ((M_(T#(i+2)))*H#(T#(i+3))) + ((M_(T#(i+4)))*H#(T#(i+5))) - ((M_(T#(i+6)))*H#(T#(i+7))))));
  P = new HashTable from apply((length G), i -> (G#i => B#i));

  return P

)

makeTableOfnxnMinors = (n, M) -> (
    rowCt := numRows M;
    colCt := numColumns M;
    gg := subsets(rowCt, n);
    ff := subsets(colCt, n);
    L := toList ((set gg)**(set ff));
    R := {};
    if (n > 2) then (R = makeTableOfnxnMinors(n-1, M)) else ( return new HashTable from  apply(L, i -> (i => det (M^(i#0)_(i#1))) ) );
    H := new HashTable from  apply(L, i -> (i => sum(0..n-1, curCol -> (-1)^(curCol)*M_(i#0#0, i#1#curCol)*R#(toList drop(i#1, {curCol, curCol}), toList drop(i#0, {0,0}))))); --det (M^(i#0)_(i#1))

    return H;
)

fastMinors = (n, M1) -> (
    --if (n == numColumns M1) then M1 = transpose M1;
    print "magic?";
    --ideal values fastMinorsTableP(n, M1)
    return 7;
);



temp = (Lis,n, M1, R1) -> (
    --print "temp started";
	K := apply(Lis, i -> (i => sum(0..n-1, curCol -> (-1)^(curCol)*M1_(i#0#(n-1), i#1#curCol)*R1#(toList drop(i#0, {n-1,n-1}), toList drop(i#1, {curCol, curCol})))));
	return K
);

fastMinorsTableP = method(Options=>{Threads=>(allowableThreads-2)});

fastMinorsTableP(ZZ, Matrix) := opts -> (n, M1) -> (
    print opts.Threads;
    if (opts.Threads <= 0) then return fastMinorsTable(n, M1);
    rowCt := numRows M1;
    colCt := numColumns M1;
    gg := subsets(rowCt, n);
    ff := subsets(colCt, n);
    L := toList ((set gg)**(set ff));
    lenL := length L;
    R1 := {};
        print "fastMinorsTableP: about to recurse";
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
    if (debugLevel > 0) then print "fastMinorsTableP: did recursion, going to do tasks";
    print "fastMinorsTableP: did recursion, going to do tasks";
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
    L := toList  ((set gg)**(set ff));
    R := {};
    if (n > 2) then (R = fastMinorsTable(n-1, submatrix'(M1, {rowCt-1}, ))) else ( return new HashTable from  apply(L, i -> (i => det (M1^(i#0)_(i#1))) ) );
    H := new HashTable from  apply(L, i -> (i => sum(0..n-1, curCol -> (-1)^(curCol)*M1_(i#0#(n-1), i#1#curCol)*R#(toList drop(i#0, {n-1,n-1}), toList drop(i#1, {curCol, curCol}))))); --det (M^(i#0)_(i#1))

    return H;
)
