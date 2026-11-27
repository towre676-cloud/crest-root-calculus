// Use a licensed Magma installation with RationalPointsGenus2.
// Retain the startup banner for the exact version. All runs use defaults.
P<x> := PolynomialRing(Rationals());
f := 6*x^5 + 28*x^3 + 54*x;
inputs := [
    <"Stoll_d26", 26*f>,
    <"Stoll_d1", f>,
    <"T1", 154*x^6-90*x^5-1605*x^4-900*x^3+5268*x^2+7200*x-1792>,
    <"T2", 28*x^6-450*x^5-1317*x^4+900*x^3+6420*x^2+1440*x-9856>,
    <"Stoll_d22", 22*f>
];
for input in inputs do
    C := HyperellipticCurve(input[2]);
    start := Cputime();
    print "Input", input[1];
    points, complete, height := RationalPointsGenus2(C);
    print "Points", points;
    print "Completeness flag", complete;
    if not complete then
        print "Certified search height", height;
    end if;
    print "CPU seconds", Cputime(start);
end for;
