#!/usr/bin/env python3
import csv,gzip,io,json,re,time,traceback,sys,os
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor,as_completed
from datetime import datetime,date
from pathlib import Path
import requests

START='2023-11-12'; END='2026-07-27'
SPECIES={
'hsa':('Homo sapiens',9606),'mmu':('Mus musculus',10090),'rno':('Rattus norvegicus',10116),
'hgl':('Heterocephalus glaber',10181),'xtr':('Xenopus tropicalis',8364),'dre':('Danio rerio',7955),
'dme':('Drosophila melanogaster',7227),'cel':('Caenorhabditis elegans',6239),
'ath':('Arabidopsis thaliana',3702),'sce':('Saccharomyces cerevisiae',4932)}
KEY=os.environ.get('SPECIES_KEY') or (sys.argv[1] if len(sys.argv)>1 else '')
if KEY not in SPECIES: raise SystemExit(f'Unknown species key: {KEY}')
SP,TAX=SPECIES[KEY]
CHUNK_TOTAL=max(1,int(os.environ.get('CHUNK_TOTAL','1')))
CHUNK_INDEX=int(os.environ.get('CHUNK_INDEX','0'))
TAG=f'{KEY}-part{CHUNK_INDEX+1}of{CHUNK_TOTAL}' if CHUNK_TOTAL>1 else KEY
OUT=Path('results')/TAG; OUT.mkdir(parents=True,exist_ok=True)
CACHE=Path('cache')/KEY; CACHE.mkdir(parents=True,exist_ok=True)
EUTIL='https://eutils.ncbi.nlm.nih.gov/entrez/eutils'; GEO='https://ftp.ncbi.nlm.nih.gov/geo/series'; ENA='https://www.ebi.ac.uk/ena/portal/api/search'
S=requests.Session(); S.headers['User-Agent']='circBank2.0-update/2.2 academic contact: circbank@example.org'
EXCL={
'polyA/oligo-dT':[r'poly\s*[-_]?\s*\(?a\)?',r'poly\s*[-_]?\s*d?t',r'oligo\s*\(?d?t\)?',r'\bmrna\s*[-_ ]?seq\b',r'mrna\s*(selection|enrichment|capture)'],
'single-cell':[r'single\s*[-_ ]?cell',r'\bscrna\b',r'single\s*[-_ ]?nucleus',r'\b10\s*x\b',r'10x\s*genomics',r'chromium',r'smart\s*[-_ ]?seq',r'cel\s*[-_ ]?seq'],
'3prime':[r"3\s*['’′-]?\s*(end|prime)",r"3\s*['’′]?\s*gene\s*expression",r'quant\s*seq',r'tag\s*[-_ ]?seq'],
'FFPE':[r'\bffpe\b',r'formalin\s*[-_ ]?fixed',r'paraffin\s*[-_ ]?embedded']}
CTRL=re.compile(r'\b(wild\s*[- ]?type|wildtype|\bwt\b|untreated|control|normal|mock|vehicle)\b',re.I)
PERT=re.compile(r'\b(knockout|knock\s*down|crispr|treated|treatment|infect|disease|tumou?r|cancer|mutant|overexpress|drug)\b',re.I)

def get(url,params=None,tries=5,timeout=180):
    err=None
    for i in range(tries):
        try:
            r=S.get(url,params=params,timeout=timeout)
            if r.status_code==200:return r
            err=RuntimeError(f'{r.status_code} {r.url} {r.text[:300]}')
        except Exception as e: err=e
        time.sleep(min(20,2**i))
    raise err

def esearch():
    term=f'"{SP}"[Organism] AND ("{START}"[PDAT] : "{END}"[PDAT]) AND "Expression profiling by high throughput sequencing"[DataSet Type]'
    r=get(f'{EUTIL}/esearch.fcgi',{'db':'gds','term':term,'retmax':100000,'retmode':'json'})
    return term,r.json()['esearchresult']['idlist']

def esummary(ids):
    out=[]
    for i in range(0,len(ids),200):
        r=get(f'{EUTIL}/esummary.fcgi',{'db':'gds','id':','.join(ids[i:i+200]),'retmode':'json'})
        js=r.json()['result']
        for uid in js.get('uids',[]):
            acc=js[uid].get('accession','')
            if acc.startswith('GSE'):out.append(acc)
        time.sleep(.34)
    return sorted(set(out))

def bucket(g):return f'GSE{int(g[3:])//1000}nnn'
def soft_url(g):return f'{GEO}/{bucket(g)}/{g}/soft/{g}_family.soft.gz'
def flat(v):return ' ; '.join(v) if isinstance(v,list) else str(v or '')

def parse_soft(g):
    cp=CACHE/f'{g}.soft.gz'
    if cp.exists(): raw=cp.read_bytes()
    else:
        raw=get(soft_url(g),timeout=300).content; cp.write_bytes(raw)
    txt=gzip.decompress(raw).decode('utf-8','replace')
    ser={'gse':g}; samples={}; cur=None
    for line in txt.splitlines():
        if line.startswith('^SERIES = '):cur=ser
        elif line.startswith('^SAMPLE = '):
            gsm=line.split('=',1)[1].strip(); samples[gsm]={'gse':g,'gsm':gsm};cur=samples[gsm]
        elif line.startswith('!') and cur is not None:
            k,_,v=line.partition(' = ');k=k[1:]
            if k in cur:
                if not isinstance(cur[k],list):cur[k]=[cur[k]]
                cur[k].append(v)
            else:cur[k]=v
    return ser,samples

def pubdate(ser):
    m=re.search(r'Public on\s+(.+)$',flat(ser.get('Series_status','')),re.I)
    if not m:return ''
    s=m.group(1).strip()
    for f in ('%b %d, %Y','%Y-%m-%d','%d-%b-%Y'):
        try:return datetime.strptime(s,f).date().isoformat()
        except:pass
    return s

def exps(sample):
    text=' ; '.join(flat(v) for v in sample.values())
    return sorted(set(re.findall(r'(?:SRX|ERX|DRX)\d+',text)))

def exclusion(text):return [lab for lab,pats in EXCL.items() if any(re.search(p,text,re.I) for p in pats)]

def ena(exp):
    fields='run_accession,experiment_accession,study_accession,secondary_study_accession,sample_accession,secondary_sample_accession,scientific_name,tax_id,first_public,library_strategy,library_source,library_selection,library_layout,instrument_platform,instrument_model,base_count,read_count,fastq_ftp,fastq_bytes'
    r=get(ENA,{'result':'read_run','query':f'experiment_accession="{exp}"','fields':fields,'format':'tsv','limit':0},timeout=180)
    return list(csv.DictReader(io.StringIO(r.text),delimiter='\t'))

def tissue(s):
    vals=s.get('Sample_characteristics_ch1',[]); vals=vals if isinstance(vals,list) else [vals]
    kv={}
    for x in vals:
        if ':' in x:
            k,v=x.split(':',1);kv[k.strip().lower()]=v.strip()
    for k in ('tissue','organ','source','cell type','cell_type','cell line','cell_line','developmental stage','stage'):
        if kv.get(k):return kv[k][:120]
    return (flat(s.get('Sample_source_name_ch1')) or flat(s.get('Sample_title')) or '未明确')[:120]

def write(name,rows):
    keys=[]
    for x in rows:
        for k in x:
            if k not in keys:keys.append(k)
    with open(OUT/name,'w',newline='',encoding='utf-8-sig') as f:
        if not keys: f.write('status\nempty\n');return
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)

def main():
    term,ids=esearch();all_gses=esummary(ids)
    gses=[g for i,g in enumerate(all_gses) if i%CHUNK_TOTAL==CHUNK_INDEX]
    (OUT/'query.json').write_text(json.dumps({'species':SP,'tax_id':TAX,'term':term,'gds_ids':len(ids),'all_gses':len(all_gses),'chunk_total':CHUNK_TOTAL,'chunk_index':CHUNK_INDEX,'chunk_gses':len(gses)},ensure_ascii=False,indent=2),encoding='utf-8')
    allser={};allsam={};soft_errors=[]
    with ThreadPoolExecutor(max_workers=10) as ex:
        fut={ex.submit(parse_soft,g):g for g in gses}
        for i,f in enumerate(as_completed(fut),1):
            g=fut[f]
            try:ser,sam=f.result();allser[g]=ser;allsam.update(sam)
            except Exception as e:soft_errors.append({'species':SP,'gse':g,'error':str(e),'url':soft_url(g)})
            if i%50==0:print(TAG,'SOFT',i,len(gses),flush=True)
    pre=[];excluded=[]
    for gsm,s in allsam.items():
        g=s['gse'];ser=allser[g];pd=pubdate(ser)
        try:d=date.fromisoformat(pd)
        except:excluded.append({'species':SP,'gse':g,'gsm':gsm,'stage':'date','reason':f'bad public date {pd}'});continue
        if not(date.fromisoformat(START)<=d<=date.fromisoformat(END)):continue
        txt=' ; '.join(flat(v) for v in list(s.values())+list(ser.values()))
        bad=exclusion(txt)
        if bad:excluded.append({'species':SP,'gse':g,'gsm':gsm,'stage':'metadata','reason':'；'.join(bad)});continue
        ee=exps(s)
        if len(ee)!=1:excluded.append({'species':SP,'gse':g,'gsm':gsm,'stage':'link','reason':f'Experiment count={len(ee)}'});continue
        pre.append((g,gsm,ee[0],ser,s,txt,pd))
    exp_list=sorted(set(x[2] for x in pre)); runmap={}
    with ThreadPoolExecutor(max_workers=20) as ex:
        fut={ex.submit(ena,x):x for x in exp_list}
        for i,f in enumerate(as_completed(fut),1):
            x=fut[f]
            try:runmap[x]=f.result()
            except Exception as e:runmap[x]=[];excluded.append({'species':SP,'srx':x,'stage':'ENA','reason':str(e)})
            if i%100==0:print(TAG,'ENA',i,len(exp_list),flush=True)
    eligible=[];review=[]
    for g,gsm,exp,ser,s,txt,pd in pre:
        runs=runmap.get(exp,[])
        if len(runs)!=1:excluded.append({'species':SP,'gse':g,'gsm':gsm,'srx':exp,'stage':'multi-run','reason':f'Run count={len(runs)}'});continue
        r=runs[0];why=[];manual=[]
        strategy=(r.get('library_strategy') or '').strip();source=(r.get('library_source') or '').strip();selection=(r.get('library_selection') or '').strip()
        if (r.get('library_layout') or '').upper()!='PAIRED':why.append('非双端')
        if (r.get('instrument_platform') or '').upper() not in ('ILLUMINA','BGISEQ'):why.append('平台不符')
        if strategy not in ('RNA-Seq','ncRNA-Seq','OTHER'):why.append('strategy不符')
        elif strategy=='OTHER':manual.append('LibraryStrategy=OTHER')
        if source.upper() not in ('TRANSCRIPTOMIC','OTHER'):why.append('source不符')
        elif source.upper()=='OTHER':manual.append('LibrarySource=OTHER')
        fq=[x for x in (r.get('fastq_ftp') or '').split(';') if x]
        if len(fq)!=2:why.append(f'FASTQ数={len(fq)}')
        try:bases=int(r.get('base_count') or 0);reads=int(r.get('read_count') or 0);rl=bases/reads if reads else 0
        except:bases=reads=0;rl=0
        if bases<=3_000_000_000:why.append('数据量≤3Gb')
        if rl<=150:why.append('平均单端读长≤150')
        base={'species':SP,'tissue':tissue(s),'gse':g,'gsm':gsm,'bioproject':r.get('secondary_study_accession') or r.get('study_accession'),'sra_study':r.get('study_accession'),'biosample':r.get('secondary_sample_accession') or r.get('sample_accession'),'srx':exp,'run':r.get('run_accession'),'geo_public_date':pd,'series_title':flat(ser.get('Series_title')),'sample_title':flat(s.get('Sample_title')),'source_name':flat(s.get('Sample_source_name_ch1')),'characteristics':flat(s.get('Sample_characteristics_ch1')),'library_strategy':strategy,'library_source':source,'library_selection':selection,'layout':r.get('library_layout'),'platform':r.get('instrument_platform'),'instrument':r.get('instrument_model'),'base_count':bases,'bases_gb':round(bases/1e9,3),'read_count':reads,'calculated_read_length':round(rl,2),'fastq_ftp':r.get('fastq_ftp'),'geo_url':f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={g}','gsm_url':f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gsm}','ena_url':f'https://www.ebi.ac.uk/ena/browser/view/{r.get("run_accession") or ""}'}
        if why:
            x=dict(base);x.update({'stage':'technical','reason':'；'.join(why)});excluded.append(x);continue
        pri='优先' if CTRL.search(txt) and not PERT.search(txt) else ('较低' if PERT.search(txt) and not CTRL.search(txt) else '中等')
        score=(bases/1e9-3)*2+(rl-150)*2+(20 if pri=='优先' else 10 if pri=='中等' else 0)+(5 if strategy!='OTHER' else -5)+(5 if source.upper()!='OTHER' else -5)
        base.update({'priority':pri,'quality_score':round(score,2),'screen_status':'人工复核' if manual else '自动通过','manual_review_reason':'；'.join(manual)})
        (review if manual else eligible).append(base)
    selected=[]; groups=defaultdict(list)
    for x in eligible:groups[(x['species'],x['tissue'])].append(x)
    for key,rows in groups.items():
        by=defaultdict(list)
        for x in rows:by[x['gse']].append(x)
        ranked=sorted(by.items(),key=lambda kv:-max(y['quality_score'] for y in kv[1]))[:3]
        for grank,(g,rr) in enumerate(ranked,1):
            for srank,x in enumerate(sorted(rr,key=lambda y:-y['quality_score'])[:3],1):
                z=dict(x);z['selection_reason']=f'{key[0]}/{key[1]}第{grank}个GSE；GSE内第{srank}个样本';selected.append(z)
    write('01_selected.csv',selected);write('02_auto_pass.csv',eligible);write('03_other_review.csv',review);write('04_excluded.csv',excluded);write('05_soft_errors.csv',soft_errors)
    summary={'species':SP,'tag':TAG,'start':START,'end':END,'GSE_discovered_all':len(all_gses),'GSE_in_chunk':len(gses),'GSE_soft_loaded':len(allser),'GEO_samples_total':len(allsam),'GEO_SRA_pre_candidates':len(pre),'auto_pass_GSE':len(set(x['gse'] for x in eligible)),'auto_pass_GSM':len(set(x['gsm'] for x in eligible)),'auto_pass_runs':len(set(x['run'] for x in eligible)),'manual_review_GSE':len(set(x['gse'] for x in review)),'manual_review_GSM':len(set(x['gsm'] for x in review)),'manual_review_runs':len(set(x['run'] for x in review)),'selected_GSE':len(set(x['gse'] for x in selected)),'selected_GSM':len(set(x['gsm'] for x in selected)),'selected_runs':len(set(x['run'] for x in selected)),'excluded_records':len(excluded),'soft_errors':len(soft_errors)}
    (OUT/'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(summary,ensure_ascii=False,indent=2))

try:main()
except Exception:
    (OUT/'fatal_error.txt').write_text(traceback.format_exc(),encoding='utf-8');print(traceback.format_exc());sys.exit(0)
