# Each corpus model through native OpenSWMM 6 (1D build); .out compared with EPA SWMM 5.2.4's (census6 e_out hash).
import json, os, subprocess, hashlib, time
B='/root/work/swmm6src/build/bin/Release/openswmm'; ROOT='/root/work/corpus'; OUT='/root/work/swmm6src/census'
recs=[json.loads(l) for l in open('/root/work/host/census6.jsonl')]
with open(f'{OUT}/census_swmm6.jsonl','w') as f:
    for r in recs:
        inp=os.path.join(ROOT,r['inp']); d=os.path.dirname(inp); n=os.path.splitext(os.path.basename(inp))[0]
        work=f'{OUT}/runs/{r["inp"].replace("/","__")}'; os.makedirs(work,exist_ok=True)
        t=time.time()
        try:
            cp=subprocess.run([B,inp,f'{work}/m.rpt',f'{work}/m.out'],capture_output=True,text=True,timeout=900,cwd=d); rc=cp.returncode
        except subprocess.TimeoutExpired: rc='timeout'
        ms=int((time.time()-t)*1000)
        h=hashlib.sha256(open(f'{work}/m.out','rb').read()).hexdigest() if os.path.exists(f'{work}/m.out') and os.path.getsize(f'{work}/m.out') else None
        rec={'inp':r['inp'],'rc':rc,'ms':ms,'out':h,'epa_out':r.get('e_out'),'epa_rc':r.get('e_rc'),'same_as_epa':h is not None and h==r.get('e_out')}
        f.write(json.dumps(rec)+'\n'); f.flush()
print('done')
