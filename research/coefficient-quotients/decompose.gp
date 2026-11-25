\\ Independent squarefree normalization (Yun in characteristic zero).
monicgcd(f,g)={my(h=gcd(f,g));h/pollead(h)};

normalize(P)={
  my(a=pollead(P), f=P/a, c=monicgcd(f, f'), w=f/c, j=1, A=1, F=a, y, z);
  while(w!=1,
    y=monicgcd(w,c); z=w/y;
    A*=z^(j\2);
    if(j%2!=0,F*=z);
    w=y; c/=y; j+=1
  );
  [A,F];
};

fixeddegree(f,d,m)={
  my(n=poldegree(f), u=x^d, q=f, r, qr, b=0, k=0, rev, coeffs);
  rev=vector(d,j,polcoef(f,n+1-j)/pollead(f));
  coeffs=Vec(Ser(Polrev(rev),x,d)^(1/m),d);
  for(j=1,d-1,u+=coeffs[j+1]*x^(d-j));
  while(q!=0,
    qr=divrem(q,u); q=qr[1]; r=qr[2];
    if(poldegree(r)>0,return([]));
    b+=r*x^k; k+=1
  );
  [u,b];
};

scan(P)={
  my(t=getabstime(), normal=normalize(P), norm=getabstime()-t, f=normal[2], n=poldegree(f), m, result, out=List());
  t=getabstime();
  for(d=2,n\3,
    if(n%d==0,
      m=n\d; result=fixeddegree(f,d,m);
      if(#result,listput(out,[d,m,Vecrev(result[1]),Vecrev(result[2])]))
    )
  );
  print("PHASES ",[norm,getabstime()-t]);
  print("SQUARE ",Vecrev(normal[1]));
  print("SOURCE ",Vecrev(f));
  for(i=1,#out,result=out[i];print("MAP ",result[1]," ",result[2],"|",result[3],"|",result[4]));
  print("DONE");
};
