from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT/'src'))
from native_windows_camb.replay import build_paper_summary

def test_replay_refuses_to_label_historical_values_as_r1_verified():
    summary=build_paper_summary({"m0":{"status":"historical_pass_needs_r1_replay","value":-395.48,"provenance":"old"}}); assert summary['publication_ready'] is False and summary['verified_results']=={} and summary['nonverified_results']['m0']['status']=='historical_pass_needs_r1_replay'
def test_replay_includes_only_verified_results_in_claim_surface():
    summary=build_paper_summary({"m0":{"status":"verified","value":-395.48,"provenance":"receipt.json"},"bench":{"status":"pending_r1","value":1.8,"provenance":"historical"}}); assert set(summary['verified_results'])=={'m0'} and set(summary['nonverified_results'])=={'bench'}
