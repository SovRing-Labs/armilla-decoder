import json, pathlib, subprocess, sys, tempfile, shutil, os
root=pathlib.Path(__file__).resolve().parents[1]; mode=sys.argv[1]
def lines(p):
 return [json.loads(l) for l in pathlib.Path(p).read_text().splitlines() if l.startswith('{')]
with tempfile.TemporaryDirectory() as t:
 t=pathlib.Path(t)
 if mode=='h264':
  r=subprocess.run([sys.executable,str(root/'media/idct_rings.py')],cwd=t,check=True,capture_output=True,text=True); print(r.stdout)
  d=json.loads((t/'idct-rings.json').read_text()); assert d['blocks']==129600 and d['clean']['mismatching_samples']==0
  assert d['fault_10pct']['mismatching_samples']>0
 elif mode=='negative':
  d=json.loads((root/'media/idct-rings.json').read_text()); assert d['ring_vs_plain_numpy_slowdown']==134.1 and d['fault_10pct']['mismatching_samples']==113 and d['fault_1pct_blocks_per_stage']['mismatching_samples']==1
 elif mode=='vp9negative':
  d=json.loads((root/'media/vp9-idct-rings/result-three-1M.json').read_text()); assert d['mismatches']==0 and d['faults_injected']==200000 and d['outputs_wrong']==17134 and d['wrong_by_group_ABC']==[0,0,17134]
 elif mode=='vp9':
  subprocess.run(['rustc','-O',str(root/'media/vp9-idct-rings/src/main.rs'),'-o',str(t/'vp9')],check=True)
  r=subprocess.run([str(t/'vp9'),'200000'],check=True,capture_output=True,text=True); print(r.stdout)
  rows=[json.loads(l) for l in r.stdout.splitlines() if l.startswith('{')]; assert rows and all(x.get('mismatches',0)==0 for x in rows)
  assert rows[0]['faults_injected']==200000 and rows[0]['outputs_wrong']>0 and sum(rows[0]['injected_by_group_ABC'])==200000
 elif mode=='jsnum':
  env=dict(os.environ,CARGO_TARGET_DIR=str(t/'target')); r=subprocess.run(['cargo','run','--offline','--release','--manifest-path',str(root/'jsnum/Cargo.toml'),'--bin','verify','--','20000'],env=env,check=True,capture_output=True,text=True); print(r.stdout)
  rows=[json.loads(l) for l in r.stdout.splitlines() if l.startswith('{')]; assert len(rows)>=10 and all(x.get('mismatches',0)==0 for x in rows)
  assert all(x.get('guard_alarms_clean',0)==0 for x in rows)
  for x in rows:
   if 'outputs_changed' in x: assert x['outputs_changed']==x['caught_by_routing_rings']
 elif mode=='fast':
  subprocess.run(['cc','-O2','-march=native',str(root/'jsnum/fast_probe.c'),'-o',str(t/'fast')],check=True)
  r=subprocess.run([str(t/'fast'),'1000000'],check=True,capture_output=True,text=True); print(r.stdout)
  d=json.loads(r.stdout); assert d['results_equal_native'] and d['guard_alarms_clean']==0 and d['faults_injected']==d['faults_caught']
 print('PASS '+mode)
