\\ Independent PARI implementation of reduction, Cantor composition and the finite witnesses.
from_coeffs(cs)=Pol(cs);
monic_bezout(a,b) = {
  my(out=gcdext(a,b),lc=pollead(out[3]));
  out/lc
};
reduced_add(left,right,F,g) = {
  my(u=left[1],v=left[2],w=right[1],z=right[2]);
  if(poldegree(u)==0,return(right)); if(poldegree(w)==0,return(left));
  my(a=monic_bezout(u,w),b=monic_bezout(a[3],v+z),d=b[3]);
  my(U=u*w/d^2,V=((b[1]*(a[1]*u*z+a[2]*w*v)+b[2]*(v*z+F))/d)%U);
  while(poldegree(U)>g,U=(F-V^2)/U;U=U/pollead(U);V=(-V)%U);
  [U,V]
};
multiple(n,D,F,g) = {
  my(R=[1,0]);while(n,if(n%2,R=reduced_add(R,D,F,g));n=n\2;if(n,D=reduced_add(D,D,F,g)));
  R
};
reduction_classes(f,p,divs) = {
  my(F=Mod(1,p)*f,g=(poldegree(f)-1)\2,odd=poldegree(f)%2,r=0,L=1,rootfound=0);
  if(!isprime(p) || p==2 || poldegree(gcd(F,deriv(F)))>0,error("bad reduction"));
  if(!odd,
    for(t=0,p-1,if(!subst(F,x,t) && sum(i=1,#divs,if(#divs[i],subst(Pol(divs[i][1]),x,Mod(t,p))==0,0))==0,
      r=Mod(t,p);rootfound=1;break));
    if(!rootfound,error("no branch chart"));L=subst(deriv(F),x,r);
    F=simplify(subst(F,x,r+L/x)*x^(2*g+2)/L^2)
  );
  if(pollead(F)!=1 || poldegree(F)!=2*g+1,error("incorrect normalized model"));
  my(points=vector(#divs,i,
    my(d=divs[i],u,v);
    if(!#d,u=Mod(1,p)*x;v=L^g,
      u=Mod(1,p)*Pol(d[1]);v=Mod(1,p)*Pol(d[2]);
      if(!odd,my(t=r+L/x);u=simplify(subst(u,x,t)*x^poldegree(u));v=simplify(subst(v,x,t)*x^(g+1)/L));
      u=u/pollead(u);v=v%u
    );
    if((F-v^2)%u,error("invalid rational divisor reduction"));
    multiple(2,[u,v],F,g)
  ));
  [F,g,points]
};
check_reductions(f,divs,ell,controls,witnesses) = {
  my(identity=[1,0]);
  for(i=1,#controls,my(tab=reduction_classes(f,controls[i][1],divs));
    for(j=1,#divs,if(multiple(controls[i][2],tab[3][j],tab[1],tab[2])!=identity,error("invalid control"))));
  for(i=1,#witnesses,my(w=witnesses[i],tab=reduction_classes(f,w[1],divs),
    images=apply(D->multiple(w[2],D,tab[1],tab[2]),tab[3]),U=0);
    for(j=1,#divs,if(images[j]!=identity,U=images[j];break));
    if(!U || multiple(ell,U,tab[1],tab[2])!=identity,error("invalid witness order"));
    for(j=1,#divs,if(multiple(w[3][j],U,tab[1],tab[2])!=images[j],error("incorrect witness equation")))
  );
  1
};
