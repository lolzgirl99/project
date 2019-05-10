fastMinors = (n, M1) -> (
    ideal values fastMinorsTableP(n, M1)
);

temp = (Lis,n, M1, R1) -> (
	K := apply(Lis, i -> (i => sum(0..n-1, curCol -> (-1)^(curCol)*M1_(i#0#(n-1), i#1#curCol)*R1#(toList drop(i#0, {n-1,n-1}), toList drop(i#1, {curCol, curCol})))));
	return K
);

fastMinorsTableP = method(Options=>{Threads=>(allowableThreads-2)});



fastMinorsTableP(ZZ, Matrix) := opts -> (n, M1) -> (
    if (opts.Threads <= 1) then (
        if (debugLevel > 0) then print "fastMinorsP: Going to single threaded version.";
        return fastMinorsTable(n, M1);
    );
    rowCt := numRows M1;
    colCt := numColumns M1;
    gg := subsets(rowCt, n);
    ff := subsets(colCt, n);
    L := toList ((set gg)**(set ff));
    lenL := length L;
    R1 := {};
    if (debugLevel > 0) then print "fastMinorsTableP: about to recurse";
    if (n > 1) then (R1 = fastMinorsTableP(n-1, submatrix'(M1, {rowCt-1}, ),opts)) else ( return entries M);
    tempL := null;
    if (debugLevel > 0) then print "fastMinorsTableP: did recursion, going to do tasks";
    taskList := apply(opts.Threads, i -> (tempL =  take(L, {(i*(lenL-1))//(opts.Threads), ((i+1)*(lenL-1))//(opts.Threads)});
        return createTask(temp, (tempL, n, M1, R1) ); ) );
    apply(taskList, t -> schedule t);
    while true do (
        if (all(taskList, t->isReady(t))) then break;
        );
    myList := flatten flatten apply(taskList, t -> taskResult(t));
    H := new HashTable from myList;
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
