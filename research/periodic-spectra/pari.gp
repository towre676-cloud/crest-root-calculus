\\ Independent PARI implementations; coefficients and operations are exact.
\\ Use PARI's default charpoly dispatcher; do not force a slow algorithm.
cr_period(a,b) = {
  my(n=#a, pr=vector(#a), k, m);
  for(i=2,n,
    k=pr[i-1];
    while(k && (a[i]!=a[k+1] || b[i]!=b[k+1]), k=pr[k]);
    if(a[i]==a[k+1] && b[i]==b[k+1], k++);
    pr[i]=k;
  );
  m=n-pr[n]; if(n%m,m=n); m;
};
cr_discriminant(a,b,m) = {
  my(p=1,q=0,r=0,s=1,np,nq,B=1,c);
  for(i=1,m,
    c=b[if(i==1,m,i-1)];
    np=(x-a[i])*p-c*r; nq=(x-a[i])*q-c*s;
    r=b[i]*p; s=b[i]*q; p=np; q=nq; B*=b[i];
  );
  [p+s,B];
};
cr_pair(D,B,n) = {
  if(n==0,return([2,D]));
  my(k=n\2, p=cr_pair(D,B,k), u=B^(2*k));
  my(odd=p[1]*p[2]-D*u);
  if(n%2,[odd,p[2]^2-2*u*B^2],[p[1]^2-2*u,odd]);
};
cr_power(D,B,n) = {
  if(n==0,return(2));
  if(n==1,return(D));
  if(n%2==0,return(cr_power(D,B,n\2)^2-2*B^n));
  my(p=cr_pair(D,B,n\2));
  p[1]*p[2]-D*B^(n-1);
};
cr_solve(a,b,method) = {
  my(L=#a,m,n,d,D,B,prev=2,cur,next,A);
  if(L<3 || #b!=L,error("bad length"));
  for(i=1,L,if(type(a[i])!="t_INT" || type(b[i])!="t_INT" || b[i]==0,error("bad coefficient")));
  if(method=="pari-dense",
    A=matrix(L,L);
    for(i=1,L,A[i,i]=a[i];n=1+i%L;A[i,n]=b[i];A[n,i]=b[i]);
    return([charpoly(A),L]);
  );
  m=cr_period(a,b); n=L/m; d=cr_discriminant(a,b,m); D=d[1];B=d[2];
  if(method=="pari-periodic-doubling",
    cur=cr_power(D,B,n),
    if(method!="pari-periodic-linear",error("unknown method"));
    cur=D;for(i=2,n,next=D*cur-B^2*prev;prev=cur;cur=next)
  );
  [cur-2*B^n,m];
};
