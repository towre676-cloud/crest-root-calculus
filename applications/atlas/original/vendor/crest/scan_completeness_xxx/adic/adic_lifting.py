"""Exact digit forests for marked cubic-quotient charts over local rings."""
import sympy as sp


class Unramified:
    """Integral polynomial basis for theta modulo a fixed monic lift."""
    def __init__(self,p,modulus,ramification=1):
        if not sp.isprime(p) or ramification not in (1,2):
            raise ValueError("Expected prime p and ramification degree one or two.")
        self.p,self.modulus,self.e=int(p),tuple(map(int,modulus)),ramification
        self.f=len(modulus)-1
        if self.f<1 or modulus[-1]!=1:
            raise ValueError("The unramified modulus must be monic of positive degree.")
        t=sp.Symbol('theta')
        if not sp.Poly.from_list(list(reversed(modulus)),t,modulus=p).is_irreducible:
            raise ValueError("The unramified modulus must be irreducible modulo p.")
        self.q=p**self.f
        self.zero=(0,)*(self.e*self.f)
        self.one=(1,)+self.zero[1:]

    def residue_lift(self,code):
        if not 0<=code<self.q:
            raise ValueError("Residue code out of range.")
        return tuple(code//self.p**i%self.p for i in range(self.f))+(0,)*((self.e-1)*self.f)

    def digit(self,value,level):
        j=level%self.e
        scale=self.p**(level//self.e)
        return sum((value[j*self.f+i]//scale%self.p)*self.p**i for i in range(self.f))

    def shift(self,value,code,level):
        vec=list(value)
        j=level%self.e
        scale=self.p**(level//self.e)
        for i in range(self.f):
            vec[j*self.f+i]+=code//self.p**i%self.p*scale
        return tuple(vec)

    def arithmetic(self,precision):
        return Arithmetic(self,precision)


class Arithmetic:
    def __init__(self,order,precision):
        self.order=order
        self.mod=order.p**((precision+order.e-1)//order.e)
        self.f,self.e=order.f,order.e
        self.zero,self.one=order.zero,order.one

    def integer(self,value):
        return (int(value)%self.mod,)+self.zero[1:]

    def add(self,a,b):
        return tuple((u+v)%self.mod for u,v in zip(a,b))
    def neg(self,a):
        return tuple(-v%self.mod for v in a)

    def mul(self,a,b):
        if a==self.zero or b==self.zero:
            return self.zero
        if a==self.one:return b
        if b==self.one:return a
        f,e,p,M=self.f,self.e,self.order.p,self.mod
        if f==1:
            if e==1:return (a[0]*b[0]%M,)
            return ((a[0]*b[0]+p*a[1]*b[1])%M,(a[0]*b[1]+a[1]*b[0])%M)
        out=[[0]*(2*f-1) for _ in range(e)]
        for i,aa in enumerate(a):
            if aa:
                for j,bb in enumerate(b):
                    if bb:
                        s=i//f+j//f
                        out[s%e][i%f+j%f]+=aa*bb*(p if s>=e else 1)
        for row in out:
            for k in reversed(range(f,len(row))):
                if row[k]:
                    for i,co in enumerate(self.order.modulus[:-1]):
                        row[k-f+i]-=co*row[k]
        return tuple(v%M for row in out for v in row[:f])

    def power(self,a,n):
        value=self.one
        while n:
            if n%2:value=self.mul(value,a)
            n//=2
            if n:a=self.mul(a,a)
        return value


class Residue:
    def __init__(self,order):
        self.order=order
        # Residues are unramified polynomials with the high pi block zero.
        self.R=order.arithmetic(1)
        self.q=order.q
        self.zero,self.one=self.R.zero,self.R.one
        self.add,self.neg,self.mul=self.R.add,self.R.neg,self.R.mul
    def inverse(self,a):
        if a==self.zero:raise ZeroDivisionError
        return self.R.power(a,self.q-2)
    def encode(self,a):
        return self.order.digit(a,0)


class System:
    """Compile rational polynomial expressions to shared arithmetic instructions."""
    def __init__(self,variables,equations):
        self.variables=tuple(variables)
        index={v:i for i,v in enumerate(variables)}
        self.ops=[]
        def visit(node):
            node=sp.sympify(node)
            if node in index:return index[node]
            if node.is_Rational: op=('constant',int(node.p),int(node.q))
            elif node.is_Add:op=('add',tuple(visit(v) for v in node.args))
            elif node.is_Mul:op=('mul',tuple(visit(v) for v in node.args))
            elif node.is_Pow and node.exp.is_Integer and node.exp>=0:op=('power',visit(node.base),int(node.exp))
            else:raise ValueError(f"Nonpolynomial term or undeclared variable: {node}")
            idx=len(variables)+len(self.ops)
            self.ops.append(op)
            index[node]=idx
            return idx
        self.outputs=[visit(e) for e in equations]
    def evaluate(self,R,point):
        if len(point)!=len(self.variables):raise ValueError("Wrong variable count.")
        vals=list(point)
        for op in self.ops:
            if op[0]=='constant':
                val=R.integer(op[1]*pow(op[2],-1,R.mod))
            elif op[0]=='power':
                val=R.power(vals[op[1]],op[2])
            else:
                ids=op[1]
                val=vals[ids[0]]
                fn=R.add if op[0]=='add' else R.mul
                for i in ids[1:]:
                    val=fn(val,vals[i])
            vals.append(val)
        return [vals[i] for i in self.outputs]


class LinearSolver:
    """Row-reduce one fixed residue Jacobian; replay row operations on new residuals."""
    def __init__(self,field,matrix,n):
        self.k,self.n=field,n
        a=[row.copy() for row in matrix]
        self.ops=[]
        self.piv=[]
        for j in range(n):
            i=len(self.piv)
            pos=next((k for k in range(i,len(a)) if a[k][j]!=field.zero),None)
            if pos is None:continue
            a[i],a[pos]=a[pos],a[i]
            self.ops.append(('swap',i,pos))
            inv=field.inverse(a[i][j])
            a[i]=[field.mul(v,inv) for v in a[i]]
            self.ops.append(('scale',i,inv))
            for r,row in enumerate(a):
                if r!=i and row[j]!=field.zero:
                    coeff=field.neg(row[j])
                    a[r]=[field.add(x,field.mul(coeff,y)) for x,y in zip(row,a[i])]
                    self.ops.append(('axpy',r,i,coeff))
            self.piv.append(j)
        self.free=[j for j in range(n) if j not in self.piv]
        self.matrix=a

    def solutions(self,rhs):
        k=self.k
        b=list(rhs)
        for op in self.ops:
            if op[0]=='swap':
                i,j=op[1:];b[i],b[j]=b[j],b[i]
            elif op[0]=='scale':
                i,coeff=op[1:];b[i]=k.mul(b[i],coeff)
            else:
                r,i,c=op[1:];b[r]=k.add(b[r],k.mul(c,b[i]))
        if any(v!=k.zero for v in b[len(self.piv):]):return
        for number in range(k.q**len(self.free)):
            digit=[k.zero]*self.n
            for j in self.free:
                number,code=divmod(number,k.q)
                digit[j]=k.order.residue_lift(code)
            for i,j in enumerate(self.piv):
                v=b[i]
                for a in self.free:
                    v=k.add(v,k.neg(k.mul(self.matrix[i][a],digit[a])))
                digit[j]=v
            yield [k.encode(v) for v in digit]


class Branch:
    """All congruence roots over one residue root, through arbitrary finite depths."""
    def __init__(self,system,order,seed):
        self.system,self.order=system,order
        self.k=Residue(order)
        self.seed=tuple(order.residue_lift(c) for c in seed)
        self.nodes=[self.seed]
        self.level=1
        R=order.arithmetic(2)
        init=system.evaluate(R,self.seed)
        if any(order.digit(v,0) for v in init):
            raise ValueError("Submitted residue seed is not a solution.")
        base=[order.digit(v,1) for v in init]
        columns=[]
        for j in range(len(seed)):
            point=list(self.seed)
            point[j]=order.shift(point[j],1,1)
            value=system.evaluate(R,point)
            columns.append([self.k.add(order.residue_lift(order.digit(v,1)),
                                       self.k.neg(order.residue_lift(b))) for b,v in zip(base,value)])
        self.linear=LinearSolver(self.k,[list(row) for row in zip(*columns)],len(seed))
    def advance(self):
        n=self.level
        R=self.order.arithmetic(n+1)
        nodes=[]
        for point in self.nodes:
            residual=self.system.evaluate(R,point)
            if any(self.order.digit(v,j) for v in residual for j in range(n)):
                raise ArithmeticError("Prior node lost a lower-order equation.")
            rhs=[self.k.neg(self.order.residue_lift(self.order.digit(v,n))) for v in residual]
            for delta in self.linear.solutions(rhs):
                nodes.append(tuple(self.order.shift(v,c,n) for v,c in zip(point,delta)))
        self.nodes=nodes
        self.level+=1


def _seeds(fa,F,marks,order):
    field=Residue(order)
    data=[]
    for hit in marks:
        P=list(map(order.residue_lift,hit['point']))
        i,j,n=hit['frame']
        inv=field.inverse(P[n])
        u,v=[field.mul(P[k],inv) for k in (i,j)]
        coords=[field.mul(q,inv) for q in P]
        def deriv(var):
            A=field.zero
            for mon,c in F.diff(fa.COORDS[var]).terms():
                value=field.R.integer(int(c.p)*pow(int(c.q),-1,order.p))
                for k,ex in enumerate(mon):value=field.mul(value,field.R.power(coords[k],ex))
                A=field.add(A,value)
            return A
        A=deriv(j)
        d=field.mul(field.neg(deriv(i)),field.inverse(A))
        g=order.residue_lift(hit['g'])
        hinv=field.inverse(field.mul(g,A))
        values={'u':field.encode(u),'v':field.encode(v),'d':field.encode(d),'g':hit['g'],
                'hinv':field.encode(hinv)}
        pp=hit['parameters']
        if hit['kind']=='I':values.update(s=pp[0],b=pp[1],c=pp[2])
        elif hit['kind']=='II':values.update(s=0,a=pp[0],b=pp[1])
        else:raise ValueError("Unknown chart in scan report.")
        values.update(zip(('e2','e1','f1','e0','f0'),hit['shape']))
        data.append((tuple(hit['frame']),hit['kind'],values))
    return data


def forest(data,scan_report,ramification=2,indices=None,include_seeds=False):
    """Yield precision, surviving class indices and total marked branch counts at each level."""
    import wide_chart as wc
    if ramification not in (1,2,"adaptive"):
        raise ValueError("Choose ramification degree 1, 2, or adaptive.")
    orders={}
    F=wc.fa.input_quartic(data)
    marks=scan_report['marks_grouped_by_quotient']
    ids=list(range(len(marks))) if indices is None else list(indices)
    if len(set(ids))!=len(ids) or any(i<0 or i>=len(marks) for i in ids):
        raise ValueError("Invalid class indices.")
    systems={}
    branches={}
    for i in ids:
        branches[i]=[]
        e=(1 if len(marks[i])==2 else 2) if ramification=="adaptive" else ramification
        if e not in orders:
            orders[e]=Unramified(scan_report['prime'],scan_report['modulus_constant_first'],e)
        order=orders[e]
        for frame,kind,values in _seeds(wc.fa,F,marks[i],order):
            key=(frame,kind)
            if key not in systems:
                variables,equations=wc.chart(F,frame,kind)
                systems[key]=System(variables,equations)
            sys=systems[key]
            branch=Branch(sys,order,[values[str(v)] for v in sys.variables])
            branch.frame,branch.kind=frame,kind
            branches[i].append(branch)
    level=1
    while True:
        counts={i:sum(len(b.nodes) for b in branches[i]) for i in ids}
        row={'precision':level,'survivors':[i for i in ids if counts[i]],'branches':counts}
        if include_seeds:
            flags=[]
            for i in ids:
                for branch in branches[i]:
                    names=list(map(str,branch.system.variables))
                    order=branch.order
                    M=order.p**((level+order.e-1)//order.e)
                    odd_mod=order.p**(level//2)
                    for node in branch.nodes:
                        coords=[node[names.index(v)] for v in ('u','v','s')]
                        if order.e==2 and any(c%odd_mod for v in coords for c in v[order.f:]):
                            continue
                        flag={'class':i,'frame':list(branch.frame),
                              'kind':'I' if branch.kind=='I' else 'IIwide','modulus':M}
                        flag.update({k:[c%M for c in v[:order.f]]
                                     for k,v in zip(('u','v','s'),coords)})
                        flags.append(flag)
            row['flags']=flags
        yield row
        for group in branches.values():
            for b in group:
                if b.nodes:b.advance()
        level+=1

if __name__=="__main__":
    import argparse
    import json
    from itertools import islice
    from pathlib import Path
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("quartic",type=Path)
    ap.add_argument("scan",type=Path,help="JSON returned by field_scan.py on the same equation")
    ap.add_argument("--depth",type=int,required=True)
    ap.add_argument("--ramification",choices=["1","2","adaptive"],default="2")
    ap.add_argument("--indices",type=int,nargs="*")
    args=ap.parse_args()
    e=args.ramification if args.ramification=="adaptive" else int(args.ramification)
    if args.depth<1:ap.error("Depth must be positive.")
    for row in islice(forest(json.loads(args.quartic.read_text()),json.loads(args.scan.read_text()),
                            ramification=e,indices=args.indices),args.depth):
        print(json.dumps(row),flush=True)
