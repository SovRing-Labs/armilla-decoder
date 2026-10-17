"""Reference-aware complex matching; compact/invariant keys never replace payload."""
from dataclasses import dataclass
import hashlib
import json
import math
import numpy as np
try:
    from .common import Refusal,ordered
except ImportError:
    from common import Refusal,ordered

OFFSETS=(-2,-1,0,1,2)
FLOOR=1e-12


@dataclass(frozen=True)
class Pilot:
    reference: str
    timebase: str
    sample_rate: int=48000
    channels: int=1
    origin: int=0
    real: float=1.
    imag: float=0.

    def validate(self):
        if type(self.reference) is not str or not self.reference or type(self.timebase) is not str or not self.timebase:raise Refusal("pilot_identity")
        if type(self.sample_rate) is not int or self.sample_rate<=0 or type(self.channels) is not int or not 1<=self.channels<=2 or type(self.origin) is not int:raise Refusal("pilot_timebase")
        if not math.isfinite(self.real) or not math.isfinite(self.imag) or abs(abs(complex(self.real,self.imag))-1)>1e-10:raise Refusal("pilot_amplitude")
        return self

    def serialize(self):
        self.validate()
        row={**vars(self),"real":float(self.real).hex(),"imag":float(self.imag).hex(),"version":1}
        return json.dumps(row,sort_keys=True,separators=(",",":")).encode("ascii")

    @classmethod
    def deserialize(cls,payload):
        try:
            d=json.loads(payload)
            if set(d)!={"reference","timebase","sample_rate","channels","origin","real","imag","version"} or d.pop("version")!=1:raise ValueError()
            d["real"]=float.fromhex(d["real"]);d["imag"]=float.fromhex(d["imag"])
            return cls(**d).validate()
        except (ValueError,TypeError,KeyError,OverflowError):raise Refusal("pilot_wire")


def _signal(x,channels):
    a=np.asarray(x,dtype=np.complex128)
    if a.ndim==1:a=a[None,:]
    if a.ndim!=2 or a.shape[0]!=channels or not 1<=a.shape[1]<=4096 or not np.all(np.isfinite(a)):raise Refusal("complex_payload_shape")
    if float(np.vdot(a,a).real)<=FLOOR:raise Refusal("silence_undefined_phase")
    return a


def _compatible(a,b):
    if not isinstance(a,Pilot) or not isinstance(b,Pilot):raise Refusal("missing_pilot")
    a.validate();b.validate()
    if a.reference!=b.reference:raise Refusal("reference_mismatch")
    if (a.timebase,a.sample_rate,a.channels)!=(b.timebase,b.sample_rate,b.channels):raise Refusal("timebase_mismatch")
    if abs(a.origin-b.origin)>2:raise Refusal("alignment_outside_budget")


def aligned_scores(a,b):
    best=None
    for lag in OFFSETS:
        shifted=np.zeros_like(b)
        if lag>=0:shifted[:,lag:]=b[:,:b.shape[1]-lag] if lag else b
        else:shifted[:,:lag]=b[:,-lag:]
        norm=float(np.linalg.norm(shifted))
        if norm<=math.sqrt(FLOOR):continue
        error=float(np.sum(abs(a/np.linalg.norm(a)-shifted/norm)**2))
        row=(error,lag)
        if best is None or row<best:best=row
    return best


def match(query,candidates,query_pilot,candidate_pilots,common_phase=False):
    """One COMMON rotation only when declared and physically supplied as a pilot.

    Alignment is bounded zero-padded sample-domain search, not circular phase
    bin reordering. Returns a proposal; exact payload remains untouched.
    """
    if type(common_phase) is not bool:raise Refusal("invariance_declaration")
    if not isinstance(query_pilot,Pilot):raise Refusal("missing_pilot")
    query_pilot.validate()
    if len(candidates)!=len(candidate_pilots) or not 1<=len(candidates)<=24:raise Refusal("candidate_budget")
    q=_signal(query,query_pilot.channels)
    original=np.asarray(query).tobytes();metadata=len(query_pilot.serialize())
    qa=q/complex(query_pilot.real,query_pilot.imag) if common_phase else q
    scores=[];lags=[]
    for candidate,pilot in zip(candidates,candidate_pilots):
        _compatible(query_pilot,pilot)
        c=_signal(candidate,pilot.channels)
        if c.shape!=q.shape:raise Refusal("candidate_shape")
        if common_phase:c=c/complex(pilot.real,pilot.imag)
        best=aligned_scores(qa,c)
        if best is None:raise Refusal("silence_undefined_phase")
        scores.append(best[0]);lags.append(best[1]);metadata+=len(pilot.serialize())
    rank=ordered(scores);minimum=scores[rank[0]]
    return {"indices":rank,"scores":scores,"lags":lags,
            "margin":None if len(rank)==1 else scores[rank[1]]-minimum,
            "ties":sum(abs(v-minimum)<=1e-12 for v in scores),
            "reference_bytes":metadata,"payload_sha256":hashlib.sha256(original).hexdigest(),
            "common_phase":bool(common_phase),"available_information":"serialized explicit pilots; exact local complex samples",
            "guarantee":"bounded normalized matching proposal; no waveform reconstruction or absolute pilot-free phase"}


def direct_complex(query,candidates):
    q=_signal(query,1)
    scores=[aligned_scores(q,_signal(c,1))[0] for c in candidates]
    return ordered(scores)


def magnitude(query,candidates):
    q=abs(np.fft.fft(query));q=q/max(np.linalg.norm(q),1e-12)
    return ordered([float(np.sum((q-abs(np.fft.fft(c))/max(np.linalg.norm(np.fft.fft(c)),1e-12))**2)) for c in candidates])


def wrapped_phase(x):
    return (float(x)+math.pi)%(2*math.pi)-math.pi
