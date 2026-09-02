from pathlib import Path
import hashlib,sys
ROOT=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT/'src'))
from native_windows_camb.integrity import sha256_file,verify_hashes

def test_sha256_file_matches_hashlib(tmp_path):
    p=tmp_path/'x.bin'; p.write_bytes(b'abc'); assert sha256_file(p)==hashlib.sha256(b'abc').hexdigest()
def test_verify_hashes_fails_closed_on_missing_file(tmp_path):
    ok,errors=verify_hashes(tmp_path,{"missing.bin":"0"*64}); assert not ok and errors==['missing: missing.bin']
