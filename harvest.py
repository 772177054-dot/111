#!/usr/bin/env python3
import csv,gzip,io,json,re,time,sys
from collections import Counter,defaultdict
from concurrent.futures import ThreadPoolExecutor,as_completed
from datetime import datetime,date
from pathlib import Path
import requests,xml.etree.ElementTree as ET

START='2023-11-12'; END='2026-07-27'
SPECIES={'Homo sapiens':9606,'Mus musculus':10090,'Rattus norvegicus':10116,'Heterocephalus glaber':10181,'Xenopus tropicalis':8364,'Danio rerio':7955,'Drosophila melanogaster':7227,'Caenorhabditis elegans':6239,'Arabidopsis thaliana':3702,'Saccharomyces cerevisiae':4932}
EUTIL='https://eutils.ncbi.nlm.nih.gov/entrez/eutils'; GEO='https://ftp.ncbi.nlm.nih.gov/geo/series'; ENA='https://www.ebi.ac.uk/ena/portal/api/search'
S=requests.Session(); S.headers['User-Agent']='circBank2.0-update/1.0 academic'
EXCL={
'polyA/oligo-dT':[r'poly\s*[-_]?\s*\(?a\)?',r'poly\s*[-_]?\s*d?t',r'oligo\s*\(?d?t\)?',r'\bmrna\s*[-_ ]?seq\b',r'mrna\s*(selection|enrichment|capture)'],
'single-cell':[r'single\s*[-_ ]?cell',r'\bscrna\b',r'single\s*[-_ ]?nucleus',r'\b10\s*x\b',r'10x\s*genomics',r'chromium',r'smart\s*[-_ ]?seq',r'cel\s*[-_ ]?seq'],
'3prime':[r"3\s*['’′-]?\s*(end|prime)",r"3\s*['’′]?\s*gene\s*expression",r'quant\s*seq',r'tag\s*[-_ ]?seq'],
'FFPE':[r'\bffpe\b',r'formalin\s*[-_ ]?fixed',r'paraffin\s*[-_ ]?embedded']}
CTRL=re.compile(r'\b(wild\s*[- ]?type|wildtype|\bwt\b|untreated|control|normal|mock|vehicle)\b',re.I)
PERT=re.compile(r'\b(knockout|knock\s*down|crispr|treated|treatment|infect|disease|tumou?r|cancer|mutant|overexpress|drug)\b',re.I)

def get(url,params=None,tries=6,timeout=180):
    err=None
    for i in range(tries):
        try:
            r=S.get(url,params=params,timeout=timeout)
            if r.status_code==200:return r
            err=RuntimeError(f'{r.status_code} {r.url} {r.text[:200]}')
        except Exception as e:err=e
        time.sleep(min(30,2**i))
    raise err

def esearch(term):
    r=get(f'{EUTIL}/esearch.fcgi',{'db':'gds','term':term,'retmax':100000,'retmode':'json'})
    return r.json()['esearchresult']['idlist']

def esummary(ids):
    out=[]
    for i in range(0,len(ids),200):
        r=get(f'{EUTIL}/esummary.fcgi',{'db':'gds','id':','.join(ids[i:i+200]),'retmode':'json'})
        js=r.json()['result']
        for uid in js.get('uids',[]):
            d=js[uid]; acc=d.get('accession','')
            if acc.startswith('GSE'):out.append(acc)
        time.sleep(.35)
    return out

def bucket(g):return f'GSE{int(g[3:])//1000}nnn'
def soft_url(g):return f'{GEO}/{bucket(g)}/{g}/soft/{g}_family.soft.gz'

def parse_soft(g):
    raw=gzip.decompress(get(soft_url(g),timeout=600).content).decode('utf-8','replace')
    ser={'gse':g}; samples={}; cur=None
    for line in raw.splitlines():
        if line.startswith('^SERIES = '):cur=ser
        elif line.startswith('^SAMPLE = '):
            gsm=line.split('=',1)[1].strip(); samples[gsm]={'gse':g,'gsm':gsm}; cur=samples[gsm]
        elif line.startswith('!') and cur is not None:
            k,_,v=line.partition(' = ');k=k[1:]
            if k in cur:
                if not isinstance(cur[k],list):cur[k]=[cur[k]]
                cur[k].append(v)
            else:cur[k]=v
    return ser,samples

def flat(v):return ' ; '.join(v) if isinstance(v,list) else str(v or '')
def pubdate(s):
    m=re.search(r'Public on\s+(.+)$',flat(s.get('Series_status','')),re.I)
    if not m:return ''
    for f in ('%b %d, %Y','%Y-%m-%d','%d-%b-%Y'):
        try:return datetime.strptime(m.group(1).strip(),f).date().isoformat()
        except:pass
    return m.group(1).strip()

def srx_from_sample(s):
    text=' ; '.join(flat(v) for v in s.values())
    return sorted(set(re.findall(r'SRX\d+',text)))

def exclusion(text):
    z=[]
    for lab,pats in EXCL.items():
        if any(re.search(p,text,re.I) for p in pats):z.append(lab)
    return z

def ena_srx(srx):
    fields='run_accession,experiment_accession,study_accession,secondary_study_accession,sample_accession,secondary_sample_accession,scientific_name,tax_id,first_public,library_strategy,library_source,library_selection,library_layout,instrument_platform,instrument_model,base_count,read_count,fastq_ftp,fastq_bytes'
    r=get(ENA,{'result':'read_run','query':f'experiment_accession="{srx}"','fields':fields,'format':'tsv','limit':0},timeout=300)
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

def main():
    out=Path('results');out.mkdir(exist_ok=True)
    gses=set(); qlog=[]
    for sp in SPECIES:
        term=f'"{sp}"[Organism] AND ("{START}"[PDAT] : "{END}"[PDAT]) AND "Expression profiling by high throughput sequencing"[DataSet Type]'
        ids=esearch(term); gs=esummary(ids);gses.update(gs);qlog.append({'species':sp,'gds_ids':len(ids),'gses':len(gs),'term':term})
        print(sp,len(ids),len(gs),flush=True)
    all_series={};all_samples={};soft_err=[]
    with ThreadPoolExecutor(max_workers=8) as ex:
        fut={ex.submit(parse_soft,g):g for g in sorted(gses)}
        for i,f in enumerate(as_completed(fut),1):
            g=fut[f]
            try:ser,sam=f.result();all_series[g]=ser;all_samples.update(sam)
            except Exception as e:soft_err.append({'gse':g,'error':str(e),'url':soft_url(g)})
            if i%100==0:print('SOFT',i,len(gses),flush=True)
    candidates=[];excluded=[]
    for gsm,s in all_samples.items():
        g=s['gse'];ser=all_series[g];pd=pubdate(ser)
        try:d=date.fromisoformat(pd)
        except:excluded.append({'gse':g,'gsm':gsm,'stage':'date','reason':f'bad public date {pd}'});continue
        if not(date.fromisoformat(START)<=d<=date.fromisoformat(END)):continue
        typ=flat(s.get('Sample_type',''))
        if typ and 'SRA' not in typ.upper():excluded.append({'gse':g,'gsm':gsm,'stage':'GEO','reason':'Sample type非SRA'});continue
        txt=' ; '.join(flat(v) for v in list(s.values())+list(ser.values()))
        bad=exclusion(txt)
        if bad:excluded.append({'gse':g,'gsm':gsm,'stage':'metadata','reason':'；'.join(bad)});continue
        srxs=srx_from_sample(s)
        if len(srxs)!=1:excluded.append({'gse':g,'gsm':gsm,'stage':'link','reason':f'SRX count={len(srxs)}'});continue
        candidates.append((g,gsm,srxs[0],ser,s,txt,pd))
    print('GEO candidates',len(candidates),flush=True)
    srxs=sorted(set(x[2] for x in candidates));runmap={}
    with ThreadPoolExecutor(max_workers=12) as ex:
        fut={ex.submit(ena_srx,x):x for x in srxs}
        for i,f in enumerate(as_completed(fut),1):
            x=fut[f]
            try:runmap[x]=f.result()
            except Exception as e:runmap[x]=[];excluded.append({'srx':x,'stage':'ENA','reason':str(e)})
            if i%200==0:print('ENA',i,len(srxs),flush=True)
    eligible=[]
    for g,gsm,srx,ser,s,txt,pd in candidates:
        runs=runmap.get(srx,[])
        if len(runs)!=1:excluded.append({'gse':g,'gsm':gsm,'srx':srx,'stage':'multi-run','reason':f'Run count={len(runs)}'});continue
        r=runs[0];why=[]
        if r.get('library_layout','').upper()!='PAIRED':why.append('非双端')
        if r.get('instrument_platform','').upper() not in ('ILLUMINA','BGISEQ'):why.append('平台不符')
        if r.get('library_strategy') not in ('RNA-Seq','ncRNA-Seq','OTHER'):why.append('strategy不符')
        if r.get('library_source','').upper() not in ('TRANSCRIPTOMIC','OTHER'):why.append('source不符')
        fq=[x for x in r.get('fastq_ftp','').split(';') if x]
        if len(fq)!=2:why.append(f'FASTQ数={len(fq)}')
        try:bases=int(r.get('base_count') or 0);reads=int(r.get('read_count') or 0);rl=bases/reads if reads else 0
        except:bases=reads=0;rl=0
        if bases<=3_000_000_000:why.append('数据量≤3Gb')
        if rl<=150:why.append('平均单端读长≤150')
        if why:excluded.append({'gse':g,'gsm':gsm,'srx':srx,'stage':'technical','reason':'；'.join(why)});continue
        pri='优先' if CTRL.search(txt) and not PERT.search(txt) else ('较低' if PERT.search(txt) and not CTRL.search(txt) else '中等')
        score=(bases/1e9-3)*2+(rl-150)*2+(20 if pri=='优先' else 10 if pri=='中等' else 0)+(5 if r.get('library_strategy')!='OTHER' else -5)
        eligible.append({'species':r.get('scientific_name'),'tissue':tissue(s),'priority':pri,'gse':g,'gsm':gsm,'bioproject':r.get('secondary_study_accession') or r.get('study_accession'),'sra_study':r.get('study_accession'),'biosample':r.get('secondary_sample_accession') or r.get('sample_accession'),'srx':srx,'run':r.get('run_accession'),'geo_public_date':pd,'series_title':flat(ser.get('Series_title')),'sample_title':flat(s.get('Sample_title')),'source_name':flat(s.get('Sample_source_name_ch1')),'characteristics':flat(s.get('Sample_characteristics_ch1')),'library_strategy':r.get('library_strategy'),'library_source':r.get('library_source'),'library_selection':r.get('library_selection'),'layout':r.get('library_layout'),'platform':r.get('instrument_platform'),'instrument':r.get('instrument_model'),'base_count':bases,'bases_gb':round(bases/1e9,3),'read_count':reads,'calculated_read_length':round(rl,2),'fastq_ftp':r.get('fastq_ftp'),'quality_score':round(score,2),'geo_url':f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={g}','gsm_url':f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gsm}','ena_url':f'https://www.ebi.ac.uk/ena/browser/view/{r.get("run_accession")}'})
    # select max 3 GSE/species+tissue and max 3 GSM/GSE
    selected=[]
    groups=defaultdict(list)
    for x in eligible:groups[(x['species'],x['tissue'])].append(x)
    for key,rows in groups.items():
        by=defaultdict(list)
        for x in rows:by[x['gse']].append(x)
        ranked=sorted(by.items(),key=lambda kv:-max(y['quality_score'] for y in kv[1]))[:3]
        for rank,(g,rr) in enumerate(ranked,1):
            for srank,x in enumerate(sorted(rr,key=lambda y:-y['quality_score'])[:3],1):
                x=dict(x);x['selection_reason']=f'{key[0]}/{key[1]}第{rank}个GSE；GSE内第{srank}个样本';selected.append(x)
    def write(name,rows):
        keys=[]
        for x in rows:
            for k in x:
                if k not in keys:keys.append(k)
        with open(out/name,'w',newline='',encoding='utf-8-sig') as f:
            w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)
    write('01_selected.csv',selected);write('02_all_eligible.csv',eligible);write('03_excluded.csv',excluded);write('04_query_log.csv',qlog);write('05_soft_errors.csv',soft_err)
    summary={'start':START,'end':END,'GSE_discovered':len(gses),'GSE_soft_loaded':len(all_series),'GEO_samples_total':len(all_samples),'GEO_SRA_candidates':len(candidates),'strict_eligible_GSE':len(set(x['gse'] for x in eligible)),'strict_eligible_GSM':len(set(x['gsm'] for x in eligible)),'strict_eligible_runs':len(set(x['run'] for x in eligible)),'selected_GSE':len(set(x['gse'] for x in selected)),'selected_GSM':len(set(x['gsm'] for x in selected)),'selected_runs':len(set(x['run'] for x in selected)),'excluded_records':len(excluded)}
    (out/'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(summary,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
