import os,sys,pathlib,json
for name in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS"):os.environ[name]="1"
sys.dont_write_bytecode=True
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parent))
from common import Refusal
checks=0
def check(condition,label):
    global checks
    checks+=1
    if not condition:raise AssertionError(label)
def refuses(fn,label,errors=(Refusal,ValueError,TypeError,OverflowError,ArithmeticError)):
    try:fn()
    except errors:return check(True,label)
    check(False,label)

import hashlib
import numpy as np
from phase import Pilot,match,direct_complex
x=np.array([1+2j,3-1j,-2+4j,5+3j,-3-1j,2+0j,4-2j,-1+1j])
p=Pilot("synthetic","block");r=Pilot("synthetic","block",real=0.,imag=1.)
before=x.tobytes()
same=match(x,[x],p,[p]);check(same["scores"][0]<1e-20,"exact complex self match")
rot=match(x,[1j*x],p,[r],common_phase=True)
check(rot["scores"][0]<1e-20,"declared supplied common phase")
check(match(x,[1j*x],p,[r],common_phase=False)["scores"][0]>.01,"phase retained unless declared")
check(rot["reference_bytes"]==len(p.serialize())+len(r.serialize()),"pilot bytes actually counted")
check(x.tobytes()==before and rot["payload_sha256"]==hashlib.sha256(before).hexdigest(),"payload exact")
check(direct_complex(x,[x,-x])[0]==0,"independent direct complex baseline")
stereo=np.stack([x,x[::-1]],axis=0);candidate=stereo*1j;candidate[1,:]*=-1
p2=Pilot("synthetic","block",channels=2);r2=Pilot("synthetic","block",channels=2,real=0.,imag=1.)
check(match(stereo,[candidate],p2,[r2],common_phase=True)["scores"][0]>.01,"relative stereo phase preserved")
temporal=x.copy();temporal[3]*=-1
check(match(x,[temporal],p,[p],common_phase=True)["scores"][0]>.001,"temporal difference retained")
refuses(lambda:match(x,[x],None,[p]),"missing pilot")
refuses(lambda:match(x,[x],p,[Pilot("synthetic","other")]),"timebase mismatch")
refuses(lambda:match(np.zeros(8,complex),[np.zeros(8,complex)],p,[p]),"silence abstains")
check(Pilot.deserialize(p.serialize())==p,"serialized pilot roundtrip")
print(json.dumps({"status":"PASS","checks":checks,"scope":"synthetic complex samples only; no acoustic model/recordings/retrieval accuracy/speed claim"}))
