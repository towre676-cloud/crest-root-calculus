// K-arithmetic for the extension from d=101 to d=200.
mb := map<PowerStructure(RngOrd)->Integers() | o :-> MinkowskiBound(o)>;
fb := map<PowerStructure(RngOrd)->Integers() | o :-> Maximum(2,BachBound(o) div 10)>;
SetClassGroupBoundMaps(fb,mb);
Q := Rationals();
K<s> := QuadraticField(-2);
for row in [[107,2,-17,4],[131,2,1,-14]] do
    a := row[2]; D := row[3]+row[4]*s;
    E := EllipticCurve([K | 3*a^2*D^2*(7+4*s),0]);
    l,h := RankBounds(E); assert h eq 0;
    G,g := TorsionSubgroup(E);
    assert #G eq 2 and g(G.1)[1] eq 0;
    print row[1],D,h,#G;
end for;
// Each final pair is the numerator and denominator of the extra abscissa.
// Denominator zero means that there is no extra finite abscissa.
for row in [[107,2,-11,10,50,1],[118,1,13,-2,1,4],[131,2,19,4,2,1],[166,1,-11,-8,0,0]] do
    a := row[2]; D := row[3]+row[4]*s;
    E := EllipticCurve([K | 3*a^2*D^2*(7+4*s),0]);
    l,h := RankBounds(E); assert l eq 1 and h eq 1;
    ok,G,mp := PseudoMordellWeilGroup(E); assert ok;
    u := map<E->ProjectiveSpace(Q,1) | [E.1,3*a*D*E.3]>;
    V,R := Chabauty(mp,u : IndexBound:=2);
    assert IsPowerOf(R,2);
    u := Extend(u); L := Codomain(u);
    expected := {L![1,0],L![0,1]};
    if row[6] ne 0 then Include(~expected,L![row[5],row[6]]); end if;
    result := {u(mp(P)) : P in V};
    assert result eq expected;
    print row[1],D,R,result;
end for;
