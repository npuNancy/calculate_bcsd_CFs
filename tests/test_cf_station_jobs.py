"""Portable full-scope Slurm generation and JSON-only completion evidence."""
from collections import Counter
import hashlib
import json
from pathlib import Path
import subprocess
from types import SimpleNamespace

import pytest

from infos.scnet_patchify_stations import create_jobs as jobs
from infos.scnet_patchify_stations import run_job as runner

SHA = 'a' * 40
CONFIG = 'b' * 64


def options(path, *extra):
    return jobs.parser().parse_args(['--jobs-dir', str(path), '--code-sha', SHA,
                                    '--config-sha256', CONFIG, *extra])


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value))
    return path


def test_full_scope_and_dependencies(tmp_path):
    manifest, scripts = jobs.build(options(tmp_path))
    rows = manifest['jobs']; extracts = [r for r in rows if r['stage'] == 'extract']
    assert len(rows) == len(scripts) == len({r['unit_id'] for r in rows}) == 3386
    assert len(extracts) == 3384
    assert set(Counter((r['model'], r['climate_scenario']) for r in extracts).values()) == {282}
    assert len(Counter((r['model'], r['climate_scenario'], r['station_scenario']) for r in extracts)) == 36
    assert sum(r['climate_scenario'] == r['station_scenario'] for r in extracts) == 1128
    assert all(r['logical_owner'] != jobs.AGGREGATOR for r in rows)
    assert rows[0]['stage'] == 'prepare' and rows[-1]['stage'] == 'publish'
    assert all(r['depends_on'] == ['station-cf-v1/prepare'] for r in extracts)
    assert rows[-1]['depends_on'] == 'all_extract_succeeded'
    assert all(r['submit_username'] is None for r in rows)


def test_all_worker_packs_identical_and_shell_valid(tmp_path):
    """Generate every worker copy, hash all files, bash-check every unique script."""
    reference = None
    workers, _ = jobs.catalog()
    for user in workers:
        directory = tmp_path/user
        jobs.main(['--jobs-dir', str(directory), '--code-sha', SHA, '--config-sha256', CONFIG])
        actual = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in directory.iterdir()}
        assert len(actual) == 3387
        if reference is None:
            reference = actual
            manifest = json.loads((directory/'manifest.json').read_text())
            for row in manifest['jobs']:
                script = directory/row['script']
                assert actual[row['script']] == row['script_sha256']
                result = subprocess.run(['bash', '-n', str(script)], capture_output=True, text=True)
                assert result.returncode == 0, result.stderr
                text = script.read_text()
                assert '/work/home/' not in text
                assert '#SBATCH --account' not in text and '#SBATCH --chdir' not in text
                assert 'source "$SCF_ENV_FILE"' in text
                assert 'SCF_JOB_SCRIPT' in text
        else:
            assert actual == reference


@pytest.mark.parametrize('extra', [('--processes', '11'), ('--cpus-per-task', '0'),
                                    ('--prepare-cpus', '0'), ('--time', '00:00:00'),
                                    ('--partition', 'x\n#SBATCH --nodes=9'),
                                    ('--resource-profile', '../bad'), ('--station-chunk', '0')])
def test_bad_resources_rejected(tmp_path, extra):
    with pytest.raises(ValueError):
        jobs.build(options(tmp_path, *extra))


def test_dry_run_and_immutable_pack(tmp_path):
    path = tmp_path/'pack'
    args = ['--jobs-dir', str(path), '--code-sha', SHA, '--config-sha256', CONFIG]
    jobs.main(args+['--dry-run'])
    assert not path.exists()
    jobs.main(args)
    before = (path/'manifest.json').read_bytes()
    with pytest.raises(SystemExit):
        jobs.main(args)
    assert (path/'manifest.json').read_bytes() == before
    with pytest.raises(SystemExit):
        jobs.main(['--jobs-dir', str(jobs.HERE/'forbidden'), '--code-sha', SHA, '--config-sha256', CONFIG])
    assert not (jobs.HERE/'forbidden').exists()
    with pytest.raises(SystemExit):
        jobs.parser().parse_args(args+['--model', 'CANESM5'])


def test_prepared_production_scope(tmp_path):
    manifest, _ = jobs.build(options(tmp_path))
    tasks = [dict(unit_id=r['unit_id'], model=r['model'], climate_scenario=r['climate_scenario'],
                  station_scenario=r['station_scenario'], tech=r['tech'], source_patch=r['patch'],
                  years=r['years']) for r in manifest['jobs'] if r['stage'] == 'extract']
    prep = dict(sample=False, output_root=str(tmp_path), implementation={'code_sha': SHA},
                capacity_semantics='snapshot_total', tasks=tasks)
    prep['identity'] = runner.digest(prep)
    path = write(tmp_path/'prepared.json', prep)
    assert len(runner.prepared_contract(path, tmp_path, SHA)['tasks']) == 3384
    prep['tasks'][-1] = prep['tasks'][0]
    prep['identity'] = runner.digest({k:v for k,v in prep.items() if k != 'identity'})
    write(path, prep)
    with pytest.raises(ValueError, match='3384'):
        runner.prepared_contract(path, tmp_path, SHA)


@pytest.fixture
def evidence_case(tmp_path):
    args = runner.parser().parse_args(['--stage','extract','--code-sha',SHA,'--config-sha256',CONFIG,
                                      '--resource-profile','station_cf_v1','--model','CANESM5',
                                      '--climate-scenario','ssp126','--station-scenario','ssp585',
                                      '--tech','wind','--patch','R02C09'])
    uid = runner.job_unit(args)
    source = {'blocks':[{'identity':f'src{i}', 'time':{'count':10}} for i in range(8)]}
    prep = {'identity':'prepared', 'tasks':[{'unit_id':uid,'station_count':2,'source_key':'source'}],
            'sources':{'source':source}}
    prepared = write(tmp_path/'prepared.json',prep)
    manifest = dict(unit_id=uid, identity='unit', prepared_identity='prepared', station_count=2,
                    implementation={'code_sha':SHA}, status='COMPLETED', blocks=[],
                    profile={'time_chunk':240,'station_chunk':1024,'compress_level':2},
                    context=dict(model=args.model,climate_scenario=args.climate_scenario,
                                 station_scenario=args.station_scenario,tech=args.tech,source_patch=args.patch))
    for i, src in enumerate(source['blocks']):
        nc = tmp_path/f'block{i}.nc'; nc.write_bytes(b'not NetCDF: JSON/stat check only')
        identity = runner.digest({'unit':'unit','source':src['identity'],'time':src['time'],'index':i})
        write(Path(str(nc)+'.json'), dict(identity=identity,status='COMPLETED',file=runner.stat(nc),
                                        source_block_identity=src['identity'],station_count=2,time_count=10))
        manifest['blocks'].append(dict(index=i,path=str(nc),identity=identity,status='COMPLETED'))
    path = write(tmp_path/'manifest.json',manifest)
    write(tmp_path/'audit.json',{'identity':'unit','status':'COMPLETED'})
    return args, {'shared':tmp_path,'prepared':prepared}, path


def test_evidence_valid_without_reading_nc(evidence_case):
    args, runtime, path = evidence_case
    value = runner.extract_evidence(path,args,runtime)
    assert value['scientific_status'] == 'COMPLETED' and len(value['artifacts']) == 8


@pytest.mark.parametrize('corruption', ['file','sidecar','audit','context','missing_block'])
def test_evidence_rejects_corruption(evidence_case, corruption):
    args, runtime, path = evidence_case
    m = runner.read(path)
    if corruption == 'file':
        Path(m['blocks'][0]['path']).write_bytes(b'changed')
    elif corruption == 'sidecar':
        p = Path(m['blocks'][0]['path']+'.json'); v=runner.read(p)
        v['source_block_identity']='other'; write(p,v)
    elif corruption == 'audit':
        write(path.parent/'audit.json', {'identity':'wrong','status':'COMPLETED'})
    elif corruption == 'context':
        m['context']['station_scenario']='ssp126'; write(path,m)
    else:
        m['blocks'].pop(); write(path,m)
    with pytest.raises(ValueError):
        runner.extract_evidence(path,args,runtime)


def test_empty_unit_requires_zero_prepared_count(evidence_case):
    args, runtime, path = evidence_case
    m=runner.read(path); m.update(status='EMPTY_NO_STATIONS',station_count=0,blocks=[])
    write(path,m); write(path.parent/'audit.json',{'status':'EMPTY_NO_STATIONS','identity':'unit'})
    with pytest.raises(ValueError):
        runner.extract_evidence(path,args,runtime)
    p=runner.read(runtime['prepared']); p['tasks'][0]['station_count']=0; write(runtime['prepared'],p)
    assert runner.extract_evidence(path,args,runtime)['scientific_status'] == 'EMPTY_NO_STATIONS'


def test_aggregator_cannot_run(monkeypatch):
    args = runner.parser().parse_args(['--stage','prepare','--code-sha',SHA,
                                      '--config-sha256',CONFIG,'--resource-profile','station_cf_v1'])
    monkeypatch.setattr(runner.pwd, 'getpwuid', lambda _:SimpleNamespace(pw_name=jobs.AGGREGATOR))
    with pytest.raises(ValueError, match='unauthorized worker'):
        runner.preflight(args)


def test_preflight_release_and_script_binding(tmp_path, monkeypatch):
    """Exercise the actual prepare preflight, then reject changed frozen files."""
    worker = jobs.catalog()[0][0]
    repo = tmp_path/'repo'; repo.mkdir()
    home = tmp_path/'homes'/jobs.AGGREGATOR; home.mkdir(parents=True)
    shared = home/'run'; shared.mkdir()
    config = write(shared/'inputs/campaign_config.json', {'models':list(jobs.MODELS)})
    config_sha = runner.sha(config)
    argv = ['--jobs-dir',str(tmp_path/'pack'),'--code-sha',SHA,'--config-sha256',config_sha]
    jobs.main(argv)
    pack = tmp_path/'pack/manifest.json'
    release = write(shared/'inputs/release.json', dict(schema='station-cf-campaign-v1',run_id='run',
                    code_sha=SHA,shared_root=str(shared),models=list(jobs.MODELS),
                    climate_scenarios=list(jobs.SCENARIOS),station_scenarios=list(jobs.SCENARIOS),
                    years='2015-2060',capacity_semantics='snapshot_total',bcsd_spot_check_accepted=True,
                    config={'path':str(config),'sha256':config_sha}))
    environment = dict(SLURM_JOB_ID='123', SLURM_CPUS_PER_TASK='10',SCF_RUN_ID='run',SCF_REPO=str(repo),
                       SCF_SHARED_ROOT=str(shared),SCF_CONFIG=str(config),SCF_RELEASE_FILE=str(release),
                       SCF_RELEASE_SHA256=runner.sha(release),SCF_PACK_MANIFEST=str(pack),
                       SCF_JOB_SCRIPT=str(tmp_path/'pack/stcf_prepare.sh'),SLURM_JOB_ACCOUNT=worker)
    for k,v in environment.items(): monkeypatch.setenv(k,v)
    monkeypatch.delenv('SLURM_ARRAY_TASK_ID',raising=False)
    monkeypatch.setattr(runner, 'ROOT', repo)
    monkeypatch.setattr(runner, 'Path', lambda value: tmp_path/'homes' if str(value)=='/work/home' else Path(value))
    monkeypatch.setattr(runner.pwd,'getpwuid',lambda _:SimpleNamespace(pw_name=worker))
    monkeypatch.setattr(runner.subprocess,'check_output',lambda command,**kwargs:
                        SHA if command[1]=='rev-parse' else '' if command[1]=='status' else 'develop-patch-grid')
    args = runner.parser().parse_args(['--stage','prepare','--code-sha',SHA,
                                      '--config-sha256',config_sha,'--resource-profile','station_cf_v1'])
    assert runner.preflight(args)['job_id'] == '123'
    script = Path(environment['SCF_JOB_SCRIPT']); original=script.read_text()
    script.write_text(original+'# changed\n')
    with pytest.raises(ValueError,match='script hash'):
        runner.preflight(args)
    script.write_text(original)
    config.write_text('{}')
    with pytest.raises(ValueError,match='configuration changed'):
        runner.preflight(args)


def test_receipt_bound_to_job_and_not_overwritten(tmp_path, monkeypatch):
    args = ['--stage','prepare','--code-sha',SHA,'--config-sha256',CONFIG,'--resource-profile','station_cf_v1']
    prepared=write(tmp_path/'prepared.json',{})
    release=write(tmp_path/'release.json',{})
    config=write(tmp_path/'config.json',{})
    # The real preflight supplies a verified hash; use actual bytes for the end-of-job check.
    args[args.index(CONFIG)] = runner.sha(config)
    runtime=dict(shared=tmp_path, unit_id='station-cf-v1/prepare', job_id='123',run_id='run',
                 user=jobs.catalog()[0][0],prepared=prepared,prepared_hash=None,release_path=release,
                 release_hash=runner.sha(release),config=config,pack_hash='pack',row={'script_sha256':'script'})
    monkeypatch.setattr(runner,'preflight',lambda _:runtime)
    monkeypatch.setattr(runner,'execute',lambda *_:{'scientific_status':'COMPLETED','units':3384})
    runner.main(args)
    receipt=tmp_path/'runtime/receipts/station-cf-v1/prepare/123.json'
    value=runner.read(receipt)
    assert value['slurm_job_id']=='123' and value['submit_username']==runtime['user']
    assert value['status']=='COMPLETED' and value['prepared_sha256']==runner.sha(prepared)
    with pytest.raises(ValueError,match='already exists'):
        runner.main(args)


def test_script_identity_with_unset_bash_source(tmp_path):
    manifest, scripts = jobs.build(options(tmp_path))
    text = scripts[manifest['jobs'][0]['script']]
    assert 'BASH_SOURCE' not in text
    assignment = next(line for line in text.splitlines() if line.startswith('SCF_JOB_SCRIPT='))
    script = tmp_path/'slurm_script'
    script.write_text('set -eu\n'+assignment+'\nprintf "%s\\n" "$SCF_JOB_SCRIPT"\n')
    result = subprocess.run(['bash','-c',script.read_text(),str(script)],capture_output=True,text=True,check=True)
    assert result.stdout.strip()==str(script.resolve())
