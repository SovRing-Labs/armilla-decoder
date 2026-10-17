"""Canonical state identity and explicit refusals for opt-in wave APIs."""
import hashlib
import math
import struct
from dataclasses import dataclass

MAX_BYTES=1_048_576


class Refusal(ValueError):
    def __init__(self, reason):
        self.reason=reason
        super().__init__(reason)


def canonical(value):
    """Typed length framing; finite binary64 bits retain -0.0; no silent coercion."""
    active=set()
    def encode(x,depth):
        if depth>32:raise Refusal("canonical_depth")
        if x is None:return b"n"
        if type(x) is bool:return b"b"+bytes([x])
        if type(x) is int:
            try:raw=str(x).encode()
            except ValueError:raise Refusal("integer_encoding_budget")
            if len(raw)>MAX_BYTES-5:raise Refusal("canonical_size")
            return b"i"+len(raw).to_bytes(4,"big")+raw
        if type(x) is float:
            if not math.isfinite(x):raise Refusal("nonfinite")
            return b"f"+struct.pack(">d",x)
        if type(x) in (str,bytes):
            if len(x)>MAX_BYTES-5:raise Refusal("canonical_size")
            raw=x.encode("utf8") if type(x) is str else x
            if len(raw)>MAX_BYTES-5:raise Refusal("canonical_size")
            return (b"s" if type(x) is str else b"y")+len(raw).to_bytes(4,"big")+raw
        if type(x) not in (tuple,list,dict):raise Refusal("unsupported_type")
        if id(x) in active:raise Refusal("cyclic_payload")
        if len(x)>MAX_BYTES-5:raise Refusal("canonical_size")
        active.add(id(x))
        try:
            if type(x) is dict:
                if any(type(k) is not str for k in x):raise Refusal("dictionary_key")
                iterator=(encode(k,depth+1)+encode(x[k],depth+1) for k in sorted(x))
                tag=b"d"
            else:
                iterator=(encode(v,depth+1) for v in x)
                tag=b"t" if type(x) is tuple else b"l"
            parts=[];total=5
            for part in iterator:
                total+=len(part)
                if total>MAX_BYTES:raise Refusal("canonical_size")
                parts.append(part)
            return tag+len(parts).to_bytes(4,"big")+b"".join(parts)
        finally:active.remove(id(x))
    result=encode(value,0)
    if len(result)>MAX_BYTES:raise Refusal("canonical_size")
    return result


def digest(value):
    return hashlib.sha256(canonical(value)).hexdigest()


@dataclass(frozen=True)
class Context:
    reference: str
    timebase: str
    source: str
    model: str
    precision: str="exact"
    rounding: str="integer-floor"
    schema: str="wave-envelope-v1"
    sample_rate: int=48000

    def token(self):
        if any(not isinstance(v,str) or not v for v in (
            self.reference,self.timebase,self.source,self.model,self.precision,self.rounding,self.schema)):
            raise Refusal("invalid_context")
        if type(self.sample_rate) is not int or self.sample_rate<=0:raise Refusal("sample_rate")
        return tuple(vars(self).items())


def ordered(scores):
    """Deterministic total order for finite scores, lower is better."""
    if any(not math.isfinite(float(s)) for s in scores):raise Refusal("nonfinite_score")
    return sorted(range(len(scores)),key=lambda i:(scores[i],i))
