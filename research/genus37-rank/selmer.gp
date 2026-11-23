\\ The caller supplies the coefficient fields, finite set S, and witnessed local points.
\\ For a place p, each component stores the prime, a uniformizer and its unit square classes.
local_data(field,p) = {
  my(primes=idealprimedec(field,p));
  vector(#primes,i,
    my(pr=primes[i], precision=if(p==2,2*pr[3]+1,1));
    [pr,nfbasistoalg(field,idealappr(field,pr)),idealstar(field,idealpow(field,pr,precision),1,2)]
  )
};
local_square_class(field,z,d) = {
  my(v=nfeltval(field,z,d[1]), unit=z/d[2]^v);
  concat([v%2],Vec(ideallog(field,unit,d[3])))
};
square_coordinates(field,z,p,data) = {
  if(!p,return(apply(s->(s<0),Vec(nfeltsign(field,z)))));
  concat(apply(d->local_square_class(field,z,d),data))
};
coordinates(pair,elems,p,data) = {
  Col(concat(vector(#pair,i,square_coordinates(fields[pair[i]],elems[i],p,if(p,data[i],[])))))
};
rational_square(z,p) = {
  if(!p,return(z>0)); if(!z,return(0));
  my(v=valuation(z,p), unit=z/p^v);
  if(v%2,return(0));
  if(p==2,lift(Mod(unit,8))==1,kronecker(lift(Mod(unit,p)),p)==1)
};

\\ A witness represents an actual divisor on y^2=f(x) over the completion.
witness(w,p,theta,pair,f) = {
  if(type(w)!="t_VEC",
    if(!rational_square(subst(f,x,w),p),error("rational x does not lift"));
    return(apply(t->w-t,theta))
  );
  if(w[1]=="branch",
    if(!(poldegree(f)%2),error("branch witness requires odd model"));
    my(h=cubics[pair[w[2]]],k=f/h);
    return(apply(t->(-1)^poldegree(h)*subst(h,x,t)+(-1)^poldegree(k)*subst(k,x,t),theta))
  );
  my(def,off);
  if(w[1]=="quadratic",
    def=if(p==2,y^2+y+1,y^2-lift(znprimroot(p)));off=1,
    if(w[1]!="extension" || p!=2,error("invalid auxiliary extension"));
    def=[y^2+1,y^2-2,y^2+2,y^2-3,y^2-6,y^2+6,y^3-2,y^3-y+1,y^4+y+1][w[2]];off=2
  );
  my(K=nfinit(def),at=idealprimedec(K,p),t=(w[off+1]+w[off+2]*Mod(y,def))/p^w[off+3]);
  if(#at!=1 || !nfislocalpower(K,at[1],subst(f,x,t),2),error("extension x does not lift"));
  my(h=charpoly(t));
  apply(a->(-1)^poldegree(h)*subst(h,x,a),theta)
};

rank2(M)=matrank(Mod(M,2));
local_rows(pair,units,p,ws) = {
  my(f=prod(i=1,#pair,cubics[pair[i]]),g=(poldegree(f)-1)\2,odd=poldegree(f)%2);
  my(data=if(p,apply(j->local_data(fields[j],p),pair),[]));
  my(theta=apply(j->Mod(x,cubics[j]),pair),one=vector(#pair,i,1));
  my(dim=#coordinates(pair,one,p,data),span=matrix(dim,0));
  if(!odd,
    my(scalars=if(!p,[-1],if(p==2,[-1,2,5],[p,lift(znprimroot(p))])));
    for(i=1,#scalars,span=concat(span,coordinates(pair,one*scalars[i],p,data)))
  );
  my(expected=rank2(span));
  if(p,expected+=sum(i=1,#pair,#data[i])-if(odd,1,2)+if(p==2,g,0),
       expected+=(sum(i=1,#pair,fields[pair[i]].sign[1])-if(odd,1,2))\2);
  for(i=1,#ws,span=concat(span,coordinates(pair,witness(ws[i],p,theta,pair,f),p,data)));
  if(rank2(span)!=expected,error("local Kummer image does not have the required rank"));
  my(U=matrix(dim,0));
  for(i=1,#pair,for(j=1,#units[i],
    U=concat(U,coordinates(pair,vector(#pair,k,if(k==i,units[i][j],1)),p,data))
  ));
  if(matsize(span)[2]==0,return(U));
  lift(matker(Mod(span~,2))~*Mod(U,2))
};

cubics=[x^3+15*x^2-34*x+1,x^3+10*x^2-19*x+4,x^3+11*x^2+3*x-14,x^3+11*x^2-11*x+16];
fields=apply(f->bnfinit(f,1),cubics);
for(i=1,4,if(!bnfcertify(fields[i]),error("uncertified class group or units")));

sector(pair,S,ws) = {
  my(f=prod(i=1,#pair,cubics[pair[i]]),odd=poldegree(f)%2);
  my(bad=select(p->valuation(poldisc(f),p)>=if(#pair==4,1,2),Vec(factor(abs(poldisc(f)))[,1])));
  if(!setisset(S) || !setsearch(S,2) || #setminus(bad,S),error("missing bad prime"));
  my(units=vector(#pair,i,
    my(K=fields[pair[i]],above=concat(apply(p->idealprimedec(K,p),S)));
    if(bnfsunit(K,above)[5][1]%2==0,error("even S-class number"));
    apply(u->nfbasistoalg(K,nffactorback(K,u)),bnfunits(K,above)[1])
  ));
  my(n=sum(i=1,#pair,#units[i]),norms=concat(apply(v->apply(z->norm(z),v),units)));
  my(R=matrix(1+#S,n,i,j,if(i==1,norms[j]<0,valuation(norms[j],S[i-1])%2)));
  for(i=1,#S+1,my(A=local_rows(pair,units,if(i==1,0,S[i-1]),ws[i]));R=matconcat([R;A]));
  my(diagonal=if(odd,0,1+#S),torsion=#pair-if(odd,1,2));
  my(bits=vector(matsize(R)[1],i,sum(j=1,n,2^(j-1)*(R[i,j]%2))));
  print([pair,n,diagonal,torsion,n-rank2(R)-diagonal-torsion,bits]);
};
