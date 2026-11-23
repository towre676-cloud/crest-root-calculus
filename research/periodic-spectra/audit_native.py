"""Replay both studies, including the after-the-fact compiled-PARI control."""
import json
from pathlib import Path

from audit import audit, audit_data, digest, require

HERE=Path(__file__).resolve().parent


def audit_native():
    original=audit()
    raw=json.loads((HERE/"native-results.json").read_text())
    protocol=json.loads((HERE/"native-results-protocol.json").read_text())
    panel=json.loads((HERE/"panel.json").read_text())
    reference=json.loads((HERE/"results-reference.json").read_text())
    require(raw["protocol_sha256"]==digest(HERE/"native-results-protocol.json"),"Supplement protocol changed")
    for key,path in [("panel_sha256","panel.json"),("reference_sha256","results-reference.json"),
                     ("original_protocol_sha256","results-protocol.json"),("original_results_sha256","results.json")]:
        require(protocol[key]==digest(HERE/path),"Supplement provenance changed: "+path)
    for name,sha in protocol["source_sha256"].items():
        require(digest(HERE/name)==sha,"Supplement source changed: "+name)
    return {"original":original,"compiled_pari_supplement":audit_data(raw,protocol,reference,panel)}


if __name__=="__main__":
    print(json.dumps(audit_native(),indent=2))
