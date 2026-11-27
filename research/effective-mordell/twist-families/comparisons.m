// For the hosted calculator, run one procedure call per submission.
// The 60-second limit on the service is separate from a false return flag.
Q<x> := PolynomialRing(Rationals());
procedure CompareCurve(f)
    t := Cputime();
    P,complete,height := RationalPointsGenus2(HyperellipticCurve(f));
    print P,complete,Cputime(t);
    if not complete then print height; end if;
end procedure;
F := 6*x^5+28*x^3+54*x;
CompareCurve(22*F);
CompareCurve(55*F);
CompareCurve(4*x^5-30*x^4-15*x^3-5*x^2-24*x);
CompareCurve(78*F);
CompareCurve(118*F);
CompareCurve(131*F);
