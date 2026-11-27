// Arithmetic over Q(sqrt(-2)). Run in Magma V2.29-10 or a compatible release.
mb := map<PowerStructure(RngOrd)->Integers() | o :-> MinkowskiBound(o)>;
fb := map<PowerStructure(RngOrd)->Integers() | o :-> Maximum(2,BachBound(o) div 10)>;
SetClassGroupBoundMaps(fb,mb);
Q := Rationals();
K<s> := QuadraticField(-2);
D := -5-2*s;
E := EllipticCurve([K | 3*D^2*(7+4*s),0]);
lo,hi := RankBounds(E);
assert lo eq 1 and hi eq 1;
ok,G,mp := PseudoMordellWeilGroup(E);
assert ok;
u := map<E -> ProjectiveSpace(Q,1) | [E.1,3*D*E.3]>;
V,R := Chabauty(mp,u : IndexBound:=2);
assert IsPowerOf(R,2);
u := Extend(u);
P1 := Codomain(u);
vals := {u(mp(P)) : P in V};
assert vals eq {P1![1,0],P1![0,1],P1![1,1]};
print lo,hi,ok,R,vals;

T := EllipticCurve([K | 945+1620*s,61290-10692*s]);
lo,hi := RankBounds(T);
assert hi eq 0;
A,mA := TorsionSubgroup(T);
assert #A eq 4;
roots := {mA(v)[1] : v in A | v ne A!0};
assert roots eq {18*s-33,-36*s+30,18*s+3};
print lo,hi,#A,roots;
