"""Diagrama de Preston da SAD usando oitavas de abundancia em log2."""
from __future__ import annotations
import csv, glob, os
import matplotlib.pyplot as plt
import numpy as np
from estilo_graficos import tons_h

ARQUIVO = "ESP_SIT.dat"
MIN_CICLOS = 2
CLASSES_PRESTON = 13

def localizar_base():
    atual=os.path.abspath(os.getcwd())
    if glob.glob(os.path.join(atual,"ciclo_*")): return atual
    candidatas=sorted(p for p in glob.glob(os.path.join(atual,"**","Simula*"),recursive=True) if os.path.isdir(p))
    if not candidatas: raise RuntimeError(f"Nenhuma pasta de simulacao em {atual}")
    return candidatas[0]

BASE=localizar_base(); SAIDA=os.path.join(BASE,"SAD_PRESTON"); os.makedirs(SAIDA,exist_ok=True)

def valor(p,prefixo): return os.path.basename(p).removeprefix(prefixo).removeprefix("=")

def configuracoes(ciclos):
    itens=set()
    for ciclo in ciclos:
        for arq in glob.glob(os.path.join(ciclo,"W*","h*","NRK=*",ARQUIVO)):
            nrk=os.path.dirname(arq); h=os.path.dirname(nrk); w=os.path.dirname(h)
            n=valor(nrk,"NRK")
            if float(n) not in (75, 100):
                itens.add((valor(w,"W"),valor(h,"h"),n))
    return sorted(itens,key=lambda x:tuple(map(float,x)))

def abundancias(arq):
    dados=np.loadtxt(arq,dtype=np.int64)
    if dados.ndim==1: dados=dados.reshape(1,-1)
    especies=dados[:,1]
    especies=especies[especies!=0]
    if not len(especies): return np.array([],dtype=np.int64)
    _,contagens=np.unique(especies,return_counts=True)
    return contagens

def coletar(W,H,n,ciclos):
    series=[]; nomes=[]
    for ciclo in ciclos:
        arq=os.path.join(ciclo,f"W{W}",f"h{H}",f"NRK={n}",ARQUIVO)
        if not os.path.exists(arq): continue
        try: contagens=abundancias(arq)
        except (OSError,ValueError,IndexError) as erro:
            print(f"Falha em {arq}: {erro}"); continue
        if len(contagens): series.append(contagens); nomes.append(os.path.basename(ciclo))
    return series,nomes

def histograma_preston(contagens):
    expoentes=np.floor(np.log2(contagens.astype(float))).astype(int)
    return np.bincount(expoentes,minlength=CLASSES_PRESTON)[:CLASSES_PRESTON]

def gerar(W,H,ns,ciclos,registros):
    grupos=[]
    for n in ns:
        series,nomes=coletar(W,H,n,ciclos)
        if len(series)>=MIN_CICLOS: grupos.append((n,series,nomes))
        else: print(f"W={W} h={H} n={n}: sem ciclo complementar")
    if not grupos: return False

    max_oitava=CLASSES_PRESTON
    fig,ax=plt.subplots(figsize=(14,7),dpi=180)
    cores=tons_h(H,len(grupos))
    shift=0.; xticks=[]; xlabels=[]; handles=[]; espaco=2.0

    for indice,((n,series,nomes),cor) in enumerate(zip(grupos,cores)):
        # Cada especie de cada ciclo entra como uma observacao no conjunto
        # agregado. O histograma de Preston e calculado uma unica vez sobre
        # esse conjunto, sem media nem desvio entre replicas.
        abundancias_agregadas=np.concatenate(series)
        contagem=histograma_preston(abundancias_agregadas).astype(float)
        pos=shift+np.arange(1,max_oitava+1)
        ax.bar(pos,contagem,width=.82,color=cor,edgecolor="#303030",linewidth=.8)
        xticks.extend(pos); xlabels.extend(map(str,range(1,CLASSES_PRESTON+1)))
        centro=shift+(max_oitava+1)/2
        ax.text(centro,-.14,rf"$n_r={n}$",transform=ax.get_xaxis_transform(),
                ha="center",va="top",fontweight="bold")
        handles.append(plt.Rectangle((0,0),1,1,facecolor=cor,edgecolor="#303030",label=rf"$n_r={n}$"))
        for oitava,numero_especies in enumerate(contagem,1):
            expoente=oitava-1
            limite_minimo=2**expoente
            limite_maximo=2**(expoente+1)-1
            registros.append([
                W,H,n,oitava,limite_minimo,limite_maximo,len(series),
                int(numero_especies),";".join(nomes),
            ])
        print(
            f"W={W} h={H} n={n}: {len(series)} ciclos agregados | "
            f"{len(abundancias_agregadas)} especies | {max_oitava} classes log2"
        )
        shift+=max_oitava
        if indice<len(grupos)-1: shift+=espaco

    ax.set_xticks(xticks); ax.set_xticklabels(xlabels,rotation=55,ha="right",fontsize=8)
    ax.set_xlim(0,shift+1); ax.set_ylim(bottom=0)
    ax.set_xlabel(r"$\log_2(n)$")
    ax.set_ylabel(r"$S(n)$")
    ax.legend(handles=handles,loc="upper center",bbox_to_anchor=(.5,-.27),
              ncol=len(handles),frameon=False)
    ax.grid(False); fig.subplots_adjust(bottom=.38,top=.96,left=.07,right=.99)
    destino=os.path.join(SAIDA,f"W{W}_h{H}_SAD_PRESTON.png")
    fig.savefig(destino,dpi=180,bbox_inches="tight"); plt.close(fig)
    print(f"Salvo: {destino}"); return True

def main():
    ciclos=sorted(glob.glob(os.path.join(BASE,"ciclo_*")))
    if not ciclos: raise RuntimeError(f"Nenhuma pasta ciclo_* em {BASE}")
    grupos={}
    for W,H,n in configuracoes(ciclos): grupos.setdefault((W,H),[]).append(n)
    registros=[]
    total=sum(gerar(W,H,ns,ciclos,registros) for (W,H),ns in grupos.items())
    csv_saida=os.path.join(SAIDA,"SAD_PRESTON_dados.csv")
    with open(csv_saida,"w",encoding="utf-8",newline="") as arq:
        w=csv.writer(arq); w.writerow([
            "W","h","n","classe_log2","N_min","N_max","ciclos_agregados",
            "numero_especies","ciclos_usados",
        ]); w.writerows(registros)
    print(f"Preston finalizado: {total} graficos | tabela: {csv_saida}"); return 0

if __name__=="__main__": raise SystemExit(main())
