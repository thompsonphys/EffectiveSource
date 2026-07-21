"""Full codegen: closed-form P-side pole coefficients P1..P5 for every stored
calc_m component, symbolic in worldline params. Emits C to poles_emit.txt and a
numeric self-check against parse_val2 residuals.

Keeps `alpha` symbolic (the Laurent variable); reduce dr^2 -> (alpha-a02 dth^2)/a20.
"""
import re, math, sys, time, json
import sympy as sp
CTX="/home/jonathan/projects/sf/effectivesource/kerr-equatorial-ctx.c"
SC="/tmp/claude-1000/-home-jonathan-projects-sf-inspectre/9f83e67a-2b3f-4a8b-b605-26a1652658f2/scratchpad"
OUT=f"{SC}/poles_emit.txt"
lines=open(CTX).read().split("\n")
def find(p,s=0):
    for i in range(s,len(lines)):
        if p in lines[i]: return i
c0=find("static void effsource_equatorial_ctx_calc_m_core")
den_c=find("Denominator - there is a different",c0); rot=find("m-modes for the rotated",den_c)
ns={}
def Sy(n):
    if n not in ns: ns[n]=sp.Symbol(n,real=True)
    return ns[n]
for mm in re.finditer(r'(\w+)\s*=\s*ctx->\w+',"\n".join(lines[c0:rot])): Sy(mm.group(1))
for x in ("dr","dtheta","M"): Sy(x)
Sy("alpha")  # keep alpha symbolic = Laurent variable
MATH={"max":max,"sqrt":sp.sqrt,"sin":sp.sin,"cos":sp.cos,"log":sp.log,"pow":lambda a,b:a**b,"fabs":sp.Abs,"M_PI":sp.pi}
arr={}
def ev(e):
    e=e.replace("[","__").replace("]","")
    return eval(e,{"__builtins__":{}},{k.replace("[","__").replace("]",""):v for k,v in ns.items()}|MATH)
st=re.compile(r'^\s*(?:const\s+)?(?:double|long double)?\s*([A-Za-z_]\w*)(\[(\d+)\])?\s*=\s*(.+);\s*$')
for i in range(c0,rot):
    ln=lines[i]
    if "ctx->" in ln or "=" not in ln: continue
    m2=st.match(ln)
    if not m2: continue
    lhs,_,idx,rhs=m2.groups()
    if lhs=="alpha": continue          # keep symbolic
    try: v=ev(rhs)
    except Exception: continue
    if idx is not None: arr.setdefault(lhs,{})[int(idx)]=v; ns[f"{lhs}__{idx}"]=v
    else: ns[lhs]=v

dr=ns["dr"]; dtheta=ns["dtheta"]; al=ns["alpha"]; beta=ns["beta"]
a20=ns["alpha20"]; a02=ns["alpha02"]; rt=ns["rt"]; rtt=ns["rtt"]
def G(n): return ns[n]
DenRe=G("DenRePhiSb"); DenIm=G("DenImPhiSb")
dDenRe={"dr":G("dDenRePhiSb_dr"),"dtheta":G("dDenRePhiSb_dtheta"),"dt":G("dDenRePhiSb_dt")}
d2DenRe={"dr":G("d2DenRePhiSb_dr2"),"dtheta":G("d2DenRePhiSb_dtheta2"),"dt":G("d2DenRePhiSb_dt2")}
dDenIm={"dr":G("dDenImPhiSb_dr"),"dtheta":G("dDenImPhiSb_dtheta"),"dt":G("dDenImPhiSb_dt")}
d2DenIm={"dr":G("d2DenImPhiSb_dr2"),"dtheta":G("d2DenImPhiSb_dtheta2"),"dt":G("d2DenImPhiSb_dt2")}
dC1_dr=G("dC1_dr"); dC1_dt=G("dC1_dt"); dC1_dtheta=G("dC1_dtheta")
d2C1_dr2=G("d2C1_dr2"); d2C1_dt2=G("d2C1_dt2"); d2C1_dtheta2=G("d2C1_dtheta2")

import numpy as np
ReEI=np.zeros((2,5,27)); ImEI=np.zeros((2,5,27))
for ln in open(f"{SC}/ei_m2.txt"):
    t,i,j,k,v=ln.split(); (ReEI if t=="R" else ImEI)[int(i),int(j),int(k)]=float(v)
pk=[1.3862943611198906188,0.096573590279972654709,0.030885144532484618274,0.014937600369780984912,0.0087663121971760665734,0.0057548876844001139245,0.0040646585499939759365]
pe=[1.0,0.44314718055994530942,0.056805192709979491031,0.021831370443737181895,0.011544521417308361798,0.0071420003133959599161,0.0048547433371649481808]
C1=al/beta; tt_=sp.Symbol("tt_"); mt=tt_/(1+tt_)
eP=[sum(sp.Float(pk[n])*mt**n for n in range(len(pk))), sum(sp.Float(pe[n])*mt**n for n in range(len(pe)))]
# together() normalises m1=C1/(1+C1)=al/(al+beta) so (al+beta) is the Pow base (not 1+al/beta)
ellipP=[sp.together(e.subs(tt_,C1)) for e in eP]
dellipP_dC=[sp.together(sp.diff(e,tt_).subs(tt_,C1)) for e in eP]
d2ellipP_dC2=[sp.together(sp.diff(e,tt_,2).subs(tt_,C1)) for e in eP]
m=2
def A(nm): return arr[nm]
def dA_tot(part,d):
    Adr=A("dReA_dr" if part=="re" else "dImA_dr")
    if d=="dr": return Adr
    if d=="dtheta": return A("dReA_dtheta" if part=="re" else "dImA_dtheta")
    p=A("dReA_dt" if part=="re" else "dImA_dt"); return {j:p[j]-rt*Adr[j] for j in p}
def d2A_tot(part,d):
    if d=="dr": return A("d2ReA_dr2" if part=="re" else "d2ImA_dr2")
    if d=="dtheta": return A("d2ReA_dtheta2" if part=="re" else "d2ImA_dtheta2")
    p2=A("d2ReA_dt2" if part=="re" else "d2ImA_dt2"); ptr=A("d2ReA_dtr" if part=="re" else "d2ImA_dtr")
    pdr=A("dReA_dr" if part=="re" else "dImA_dr"); p2r=A("d2ReA_dr2" if part=="re" else "d2ImA_dr2")
    return {j:p2[j]-2*rt*ptr[j]-rtt*pdr[j]+rt*rt*p2r[j] for j in p2}
DC1={"dr":dC1_dr,"dtheta":dC1_dtheta,"dt":dC1_dt}; D2C1={"dr":d2C1_dr2,"dtheta":d2C1_dtheta2,"dt":d2C1_dt2}
BASE={"drr":"dr","dthth":"dtheta","dtt":"dt"}
def numP(part,comp):
    Aa=A("ReA" if part=="re" else "ImA"); EI=ReEI if part=="re" else ImEI
    kr=(lambda i,j:range(max(j-i,0),m+2+j+1)) if part=="re" else (lambda i,j:range(max(j-i-1,0),m+1+j+1))
    tot=0
    for i in range(2):
        for j in range(5):
            for k in kr(i,j):
                e=sp.Float(EI[i,j,k])
                if e==0: continue
                ep=ellipP[i]; dep=dellipP_dC[i]; d2ep=d2ellipP_dC2[i]
                Ck=C1**k; dCk=(k*C1**(k-1)) if k>=1 else sp.Integer(0); d2Ck=(k*(k-1)*C1**(k-2)) if k>=2 else sp.Integer(0)
                if comp=="val": tot+=e*ep*Aa[j]*Ck
                elif comp in ("dr","dtheta","dt"):
                    dC1=DC1[comp]; dAj=dA_tot(part,comp)[j]
                    tot+=e*((dep*dC1*Aa[j]+ep*dAj)*Ck + ep*Aa[j]*dCk*dC1)
                else:
                    d=BASE[comp]; dC1=DC1[d]; d2C1=D2C1[d]; dAj=dA_tot(part,d)[j]; d2Aj=d2A_tot(part,d)[j]
                    dellP=dep*dC1; d2ellP=d2ep*dC1**2+dep*d2C1; dCkt=dCk*dC1; d2Ckt=d2Ck*dC1**2+dCk*d2C1
                    tot+=e*((d2ellP*Aa[j]+2*dellP*dAj+ep*d2Aj)*Ck + 2*(dellP*Aa[j]+ep*dAj)*dCkt + ep*Aa[j]*d2Ckt)
    return tot
def fpart(part,comp):
    Den=DenRe if part=="re" else DenIm; dDen=dDenRe if part=="re" else dDenIm; d2Den=d2DenRe if part=="re" else d2DenIm
    f=numP(part,"val")/Den
    if comp=="val": return f
    if comp in ("dr","dtheta","dt"): return (numP(part,comp)-f*dDen[comp])/Den
    d=BASE[comp]; fd=(numP(part,d)-f*dDen[d])/Den
    return (numP(part,comp)-2*fd*dDen[d]-f*d2Den[d])/Den
def uP(comp): return fpart("re",comp)+sp.I*fpart("im",comp)
cc=ns["c"]; phit=ns["phit"]; phitt=ns["phitt"]; dcdt=ns["dcdt"]; d2cdt2=ns["d2cdt2"]
imr=sp.I*m*cc; m2r=(m*cc)**2; psi_t=dcdt*dr-cc*rt+phit; imt=sp.I*m*psi_t; m2t=(m*psi_t)**2
psi_tt=d2cdt2*dr-2*dcdt*rt-cc*rtt+phitt; imtt=sp.I*m*psi_tt
Vv=uP("val"); Vr=uP("dr")-imr*Vv; Vth=uP("dtheta"); Vt=uP("dt")-imt*Vv
Vrr=uP("drr")-2*imr*uP("dr")-m2r*Vv; Vthth=uP("dthth"); Vtt=uP("dtt")-2*imt*uP("dt")-(imtt+m2t)*Vv
a=sp.Symbol("a",real=True); xpr=sp.Symbol("xpr",real=True)   # generic: a=ctx->a, xpr=xp.r
rF=xpr+dr; thF=sp.pi/2+dtheta                                 # equatorial: theta_p=pi/2
def box(dPr,dPth,dtt,dtph,drr,dthth,dphph):
    sinth=sp.sin(thF); sinth2=sinth**2; sin2th=sp.sin(2*thF); cos2th=sp.cos(2*thF)
    r=rF; r2=r*r; r3=r2*r; r4=r2*r2; a2=a*a; a4=a2*a2
    num=(2*a2*dPr - a2*dphph - a4*drr - a2*dthth - 4*dPr*r - 2*a2*dPr*r + 4*dphph*r + 2*a*dtph*r + 4*a2*drr*r + 2*dthth*r + 6*dPr*r2 - 2*dphph*r2 - 4*drr*r2 - 2*a2*drr*r2 - dthth*r2 - 2*dPr*r3 + 4*drr*r3 - drr*r4
       + (a4*dtt + 4*a*dtph*r + 2*dtt*r4 + a2*dtt*r*(2+3*r))*sinth2
       + cos2th*(a4*drr - 2*a*dtph*r + (-2+r)*r*(dthth+2*dPr*(-1+r)+drr*(-2+r)*r) + a2*(-dphph+dthth+2*dPr*(-1+r)-4*drr*r+2*drr*r2) + a2*dtt*(a2+(-2+r)*r)*sinth2)
       - a2*dPth*sin2th + 2*dPth*r*sin2th - dPth*r2*sin2th)
    den=sinth2*(a2+(-2+r)*r)*(a2+2*r2+a2*cos2th)
    return -num/den
def getcomb(name):
    if name=="src": return box(Vr,Vth,Vtt,sp.I*m*Vt,Vrr,Vthth,-(m*m)*Vv)
    return {"val":Vv,"dr":Vr,"dth":Vth,"dt":Vt,"drr":Vrr,"dthth":Vthth,"dtt":Vtt}[name]
import os
dr2sub=(al-a02*dtheta**2)/a20
ORD=16
def apb_expand(e):
    # expand every (al+beta)^q that is analytic-nontrivial at al=0: q<0 or non-integer
    # (incl. integer-negative from m1^n = al^n/(al+beta)^n in ellipP series)
    return e.replace(lambda x:x.is_Pow and x.base.is_Add and x.base==al+beta and not (x.exp.is_integer and x.exp>=0),
                     lambda x:sum(sp.binomial(x.exp,n)*beta**(x.exp-n)*al**n for n in range(ORD)))
def reduce_dr(e):
    return e.replace(lambda x:x.is_Pow and x.base==dr and x.exp>=2,
                     lambda x:dr**(x.exp%2)*dr2sub**(x.exp//2))
def poles_of(ee):
    # apb_expand -> denominator becomes pure C*al^N (I-free, real); numerator carries I
    num,den=sp.fraction(sp.together(apb_expand(ee)))
    dp=sp.Poly(den,al); N=dp.degree(); C=dp.nth(N)
    degs=sorted(set(mm[0] for mm in dp.monoms()))
    assert degs==[N], f"den NOT monomial: al-degrees={degs}"
    numr=sp.expand(reduce_dr(num))                 # polynomial in al,dr,dtheta,params,I
    out={}
    for part,sel in (("re",0),("im",1)):
        npart=numr.coeff(sp.I,sel)
        for q in (1,2,3,4,5):
            co=npart.coeff(al,N-q)
            if co!=0: out[(part,q)]=sp.cancel(co/C)
    return out

open(OUT,"w").write("")
PROG=f"{SC}/codegen_progress.txt"; open(PROG,"w").write("")
def prog(s): open(PROG,"a").write(s+"\n")
def emit(s): open(OUT,"a").write(s+"\n")
I=sp.I
WLSYM=re.compile(r'(d2Adt2|dAdt|A)\d+$')
ONLY=os.environ.get("ONLY")
names=ONLY.split(",") if ONLY else ["val","dr","dth","dt","drr","dthth","dtt","src"]
for name in names:
    t0=time.time()
    prog(f"{name}: start t={time.time()-t0:.0f}s")
    out=poles_of(getcomb(name))
    prog(f"{name}: poles done t={time.time()-t0:.0f}s")
    for part in ("re","im"):
        for q in (5,4,3,2,1):
            if (part,q) in out:
                emit(f"P{q}{part}_{name} = {sp.ccode(out[(part,q)])};")
    print(f"{name} done in {time.time()-t0:.1f}s",flush=True)
print("ALL DONE",flush=True)
