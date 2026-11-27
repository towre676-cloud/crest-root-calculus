Q:=Rationals(); F<s>:=QuadraticField(-3); om:=(-1+s)/2;
mb:=map<PowerStructure(RngOrd)->Integers()|o:->MinkowskiBound(o)>;
fb:=map<PowerStructure(RngOrd)->Integers()|o:->Maximum(2,BachBound(o) div 10)>;
SetClassGroupBoundMaps(fb,mb);
D:=[19,95,182]; kk:=[1,15,2]; tt:=[Q|1/4,1,4]; P<T>:=PolynomialRing(F);
for j in [1..3] do
  if j le 2 then
    h:=[F|3+5*om,120+15*om][j]; c:=-s*h/kk[j]; f:=T*(T+3*h)*(T-h);
  else
    deltaC:=4+6*om; eta:=6+8*om; c:=deltaC*s*eta;
    f:=T*(T-2*c*s)*(T+2*c/s);
  end if;
  assert IsSquare(Evaluate(f,c*tt[j]));
  E:=EllipticCurve(f); l,u:=RankBounds(E);
  ok,G,m:=PseudoMordellWeilGroup(E);
  assert l eq 1 and u eq 1 and ok;
  g:=map<E->ProjectiveSpace(Q,1)|[E.1,c*E.3]>;
  V,R:=Chabauty(m,g : IndexBound:=2);
  assert IsPowerOf(R,2);
  g:=Extend(g); L:=Codomain(g);
  vals:={g(m(z)):z in V};
  expected:=[{L![1,0],L![0,1],L![4,1],L![1,4]},
             {L![1,0],L![0,1],L![1,1],L![225,1]},
             {L![1,0],L![0,1],L![4,1],L![-1,1]}][j];
  assert vals eq expected;
  print D[j],l,u,ok,R,vals;
end for;
