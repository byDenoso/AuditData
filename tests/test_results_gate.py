from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT/'src'))
from native_windows_camb.results import publication_ready

def test_publication_ready_requires_all_required_gates_verified():
    ready,missing=publication_ready([{"gate_id":"G1","required_for_paper":True,"status":"verified"},{"gate_id":"G2","required_for_paper":True,"status":"pending_r1"}]); assert not ready and missing==['G2']
def test_publication_ready_ignores_optional_gate():
    ready,missing=publication_ready([{"gate_id":"G1","required_for_paper":True,"status":"verified"},{"gate_id":"Gx","required_for_paper":False,"status":"pending_r1"}]); assert ready and missing==[]
