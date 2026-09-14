"""Whittaker em barras: varios NRK no mesmo eixo com shift."""
from __future__ import annotations
import glob, os
import matplotlib.pyplot as plt
import numpy as np
from estilo_graficos import curva_suave, tons_h

ARQUIVO = "ESP_SIT.dat"
MIN_CICLOS = 2

def base_simulacao():
    atual = os.path.abspath(os.getcwd())
    if glob.glob(os.path.join(atual, "ciclo_*")): return atual
    candidatas = sorted(p for p in glob.glob(os.path.join(atual, "**", "Simula*"), recursive=True) if os.path.isdir(p))
    if not candidatas: raise RuntimeError(f"Nenhuma pasta de simulacao em {atual}")
    return candidatas[0]

BASE = base_simulacao()
SAIDA = os.path.join(BASE, "Whittaker")
os.makedirs(SAIDA, exist_ok=True)

def valor(p, prefixo): return os.path.basename(p).removeprefix(prefixo).removeprefix("=")

def configuracoes(ciclos):
    itens = set()
    for ciclo in ciclos:
        for arq in glob.glob(os.path.join(ciclo, "W*", "h*", "NRK=*", ARQUIVO)):
            nrk, h, w = os.path.dirname(arq), os.path.dirname(os.path.dirname(arq)), os.path.dirname(os.path.dirname(os.path.dirname(arq)))
            n = valor(nrk,"NRK")
            if float(n) not in (75, 100):
                itens.add((valor(w,"W"), valor(h,"h"), n))
    return sorted(itens, key=lambda x: tuple(map(float, x)))

def ler_rank(arq):
    dados = np.loadtxt(arq, dtype=np.int64)
    if dados.ndim == 1: dados = dados.reshape(1, -1)
    especies = dados[:, 1]
    especies = especies[especies != 0]
    if not len(especies): return np.array([])
    _, n = np.unique(especies, return_counts=True)
    return np.sort(n.astype(float) * 100 / n.sum())[::-1]

def coletar(W, H, NRK, ciclos):
    series, nomes = [], []
    for ciclo in ciclos:
        arq = os.path.join(ciclo, f"W{W}", f"h{H}", f"NRK={NRK}", ARQUIVO)
        if not os.path.exists(arq): continue
        try: r = ler_rank(arq)
        except (OSError, ValueError, IndexError) as erro:
            print(f"Falha em {arq}: {erro}"); continue
        if len(r): series.append(r); nomes.append(os.path.basename(ciclo))
    return series, nomes

def media_ranks(series):
    matriz = np.full((len(series), max(map(len, series))), np.nan)
    for i, s in enumerate(series): matriz[i, :len(s)] = s
    n = np.sum(~np.isnan(matriz), axis=0)
    validos = np.flatnonzero(n >= MIN_CICLOS)
    if not len(validos): return None
    fim = validos[-1] + 1
    return np.nanmean(matriz[:, :fim], axis=0), np.nanstd(matriz[:, :fim], axis=0, ddof=1), n[:fim]

def ticks_rank(n):
    return np.unique(np.clip([1, n//4, n//2, 3*n//4, n], 1, n))

def gerar(W, H, nrks, ciclos):
    grupos = []
    for nrk in nrks:
        series, nomes = coletar(W,H,nrk,ciclos)
        if len(series) < MIN_CICLOS:
            print(f"W={W} h={H} NRK={nrk}: sem ciclo complementar"); continue
        resumo = media_ranks(series)
        if resumo is not None: grupos.append((nrk,nomes,resumo))
    if not grupos: return False
    fig, ax = plt.subplots(figsize=(11,7), dpi=180)
    estilos = ["-", "--", "-.", ":", "-", "--"]
    cores = tons_h(H, len(grupos))
    maior_rank = 0
    for k, ((nrk,nomes,(media,desvio,nrep)),cor) in enumerate(zip(grupos,cores)):
        ranks=np.arange(1,len(media)+1); maior_rank=max(maior_rank,len(media))
        xs,ys=curva_suave(ranks,media)
        ax.plot(xs,ys,color=cor,linestyle=estilos[k % len(estilos)],
                linewidth=1.7,label=rf"$n_r={nrk}$ ({len(nomes)} ciclos)")
        ax.plot(ranks,media,"o",color=cor,markersize=2.5)
        print(f"W={W} h={H} NRK={nrk}: {len(media)} ranks; replicas/rank={nrep.min()}-{nrep.max()}")
    ax.set_xlim(1,maior_rank); ax.set_ylim(bottom=0)
    ax.set_xlabel("Rank das especies"); ax.set_ylabel("Abundancia relativa media (%)")
    ax.set_title(f"Diagrama de Whittaker — W={W} | H={H}",pad=18,fontsize=15)
    ax.legend(loc="upper center",bbox_to_anchor=(.5,-.13),
              ncol=max(1,len(grupos)),frameon=False)
    ax.grid(False); fig.subplots_adjust(bottom=.20)
    destino=os.path.join(SAIDA,f"W{W}_h{H}_WHITTAKER_CLASSICO.png")
    fig.savefig(destino,bbox_inches="tight",dpi=180); plt.close(fig); print(f"Salvo: {destino}")
    return True

def main():
    ciclos=sorted(glob.glob(os.path.join(BASE,"ciclo_*")))
    if not ciclos: raise RuntimeError(f"Nenhuma pasta ciclo_* em {BASE}")
    grupos={}
    for W,H,NRK in configuracoes(ciclos): grupos.setdefault((W,H),[]).append(NRK)
    n=sum(gerar(W,H,nrks,ciclos) for (W,H),nrks in grupos.items())
    print(f"Whittaker finalizado: {n} grafico(s) classico(s)."); return 0

if __name__ == "__main__": raise SystemExit(main())
