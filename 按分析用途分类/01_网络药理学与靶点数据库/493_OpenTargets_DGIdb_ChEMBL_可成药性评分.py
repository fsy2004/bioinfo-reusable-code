# =============================================================================
# 493_OpenTargets_DGIdb_ChEMBL_可成药性评分.py
# 用途    : 给一组靶基因打"可成药性"综合分 (临床期 / tractability / 已知互作 / 机制证据)
# 来源    : OpenTargets Platform GraphQL  https://api.platform.opentargets.org/api/v4/graphql
#           DGIdb 5.0 GraphQL             https://dgidb.org/api/graphql
#           ChEMBL REST                   https://www.ebi.ac.uk/chembl/api/data
# 补库依据 : 覆盖矩阵 cat01 仅 CTD/Swiss/GeneCards/OMIM；论文1 (IJMS 2026, 离子通道CRC)
#           用 OpenTargets+DGIdb+ChEMBL 的加权 DrugEvidenceScore 选可成药 hub
#           → 直接满足干法共病打法的 "actionable target" 硬指标。
# 依赖    : pip install requests pandas mygene        # 勿自动装 —— 先确认
# 输入    : genes = ["KCNQ2","RIPK2",...]  (HGNC symbol)
# 输出    : DataFrame[symbol, max_clinical_phase, tractability, n_dgidb, n_chembl_mech, DrugEvidenceScore]
# =============================================================================
import requests, pandas as pd, time, mygene

OT     = "https://api.platform.opentargets.org/api/v4/graphql"
DGIDB  = "https://dgidb.org/api/graphql"
CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"

def sym2ensembl(genes):
    out = mygene.MyGeneInfo().querymany(genes, scopes="symbol",
            fields="ensembl.gene", species="human", verbose=False)
    d = {}
    for o in out:
        eg = o.get("ensembl");  eg = eg[0] if isinstance(eg, list) else eg
        if eg: d[o["query"]] = eg["gene"]
    return d

def opentargets(ensg):
    q = """query($id:String!){ target(ensemblId:$id){
            tractability{ value } knownDrugs{ rows{ phase } } } }"""
    t = (requests.post(OT, json={"query":q,"variables":{"id":ensg}}, timeout=30)
            .json().get("data") or {}).get("target") or {}
    phases = [r["phase"] for r in (t.get("knownDrugs") or {}).get("rows",[]) if r.get("phase") is not None]
    tract  = 1 if any(x.get("value") for x in (t.get("tractability") or [])) else 0
    return (max(phases) if phases else 0), tract

def dgidb(genes):
    q = """query($n:[String!]){ genes(names:$n){ nodes{ name
            interactions{ drug{ name } interactionScore } } } }"""
    nodes = (((requests.post(DGIDB, json={"query":q,"variables":{"n":genes}}, timeout=30)
              .json().get("data") or {}).get("genes") or {}).get("nodes") or [])
    return {n["name"]: len(n.get("interactions",[])) for n in nodes}

def chembl_mech(symbol):
    try:
        tr = requests.get(f"{CHEMBL}/target/search?q={symbol}&format=json", timeout=30).json()
        n = 0
        for t in tr.get("targets",[])[:1]:
            m = requests.get(f"{CHEMBL}/mechanism?target_chembl_id={t['target_chembl_id']}&format=json", timeout=30).json()
            n += m.get("page_meta",{}).get("total_count",0)
        return n
    except Exception:
        return 0

def score(genes):
    ens, dg, rows = sym2ensembl(genes), dgidb(genes), []
    for g in genes:
        phase, tract = opentargets(ens[g]) if g in ens else (0,0)
        nd, nc = dg.get(g,0), chembl_mech(g); time.sleep(0.2)
        S = 0.50*(phase/4) + 0.25*tract + 0.15*min(nd,5)/5 + 0.10*min(nc,5)/5   # 论文权重
        rows.append(dict(symbol=g, max_clinical_phase=phase, tractability=tract,
                         n_dgidb=nd, n_chembl_mech=nc, DrugEvidenceScore=round(S,3)))
    return pd.DataFrame(rows).sort_values("DrugEvidenceScore", ascending=False)

if __name__ == "__main__":
    genes = ["KCNQ2","RIPK2","GALK1","LSM7"]      # TODO: 换成你的 hub 基因
    df = score(genes); print(df.to_string(index=False))
    # df.to_csv("druggability_scores.csv", index=False)
