"""Publication consumes verified receipts without touching CF output files."""
from pathlib import Path
from types import SimpleNamespace
import json

import pytest

from station_cf_fixtures import build
from prepare_station_cf import prepare
from extract_station_cf import extract_unit
from grid_cf_io import atomic_json,read_json
from station_cf_catalog import sha
from infos.scnet_patchify_stations.run_job import extract_evidence
from station_cf_publish import publish


def case(tmp_path):
    cfg=build(tmp_path,blocks=8,stations=('ssp585',))
    prepared=prepare(cfg,tmp_path/'out',sample=True);prep=read_json(prepared);rows={}
    for i,t in enumerate(prep['tasks']):
        args=SimpleNamespace(stage='extract',model=t['model'],climate_scenario=t['climate_scenario'],
             station_scenario=t['station_scenario'],tech=t['tech'],patch=t['source_patch'],
             code_sha=prep['implementation']['code_sha'],time_chunk=240,station_chunk=1024,compress_level=2)
        path=extract_unit(prepared,model=args.model,climate_scenario=args.climate_scenario,
             station_scenario=args.station_scenario,tech=args.tech,patch=args.patch,processes=1)
        rc=extract_evidence(path,args,dict(shared=tmp_path/'out',prepared=prepared))
        fixed=dict(run_id='test',unit_id=t['unit_id'],stage='extract',submit_username='worker',
             slurm_job_id=str(i+1),code_sha=args.code_sha,config_sha256='config',release_sha256='release',
             prepared_sha256=sha(prepared),resource_profile='v1',pack_sha256='pack',script_sha256='script')
        receipt=tmp_path/f'receipt{i}.json';atomic_json(receipt,dict(rc,**fixed,status='COMPLETED'))
        rows[t['unit_id']]=dict(fixed,status='succeeded',verified_at='checked',receipt_path=str(receipt),
             receipt_sha256=sha(receipt),expected_manifest=str(path.relative_to(tmp_path/'out')),
             accounting=[dict(job_id=j,state='COMPLETED',exit='0:0') for j in (str(i+1),str(i+1)+'.batch')])
    ledger=tmp_path/'ledger.json';atomic_json(ledger,dict(units=rows))
    return prepared,ledger


def test_receipt_publication_does_not_access_outputs(tmp_path,monkeypatch):
    prepared,ledger=case(tmp_path)
    original_open=Path.open;original_stat=Path.stat
    def check(path):
        if str(tmp_path/'out/outputs') in str(path):
            raise AssertionError('publication touched a CF output')
    def guarded_open(path,*a,**kw):check(path);return original_open(path,*a,**kw)
    def guarded_stat(path,*a,**kw):check(path);return original_stat(path,*a,**kw)
    monkeypatch.setattr(Path,'open',guarded_open);monkeypatch.setattr(Path,'stat',guarded_stat)
    index=read_json(publish(prepared,ledger,processes=2))
    assert index['summary']['year_nc_files']==16
    assert index['summary']['output_files_reopened']==0
    assert index['summary']['verification']=='previously_verified_extraction_receipts'
    assert len(index['entries'])==2 and all(len(e['blocks'])==8 for e in index['entries'])


@pytest.mark.parametrize('corruption',['unverified','slurm','receipt','scope','prepared'])
def test_receipt_publication_rejects_invalid_evidence(tmp_path,corruption):
    prepared,ledger=case(tmp_path);v=read_json(ledger);row=next(iter(v['units'].values()))
    if corruption=='unverified':row['status']='active'
    elif corruption=='slurm':row['accounting']=[]
    elif corruption=='receipt':Path(row['receipt_path']).write_text('{}')
    elif corruption=='scope':v['units'].pop(row['unit_id'])
    else:row['prepared_sha256']='changed'
    atomic_json(ledger,v)
    with pytest.raises(ValueError):publish(prepared,ledger,processes=2)
    assert not (tmp_path/'out/index/authoritative_index.json').exists()
