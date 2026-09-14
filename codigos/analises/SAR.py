import glob
import math
import os

import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
from estilo_graficos import cor_n, curva_suave


ARQUIVO_ALVO = "ESP_SIT.dat"
L_FIXO = 512
IGNORAR_ZERO = True
LS = [2, 4, 8, 16, 32, 64, 128, 256,512]
N_SAMPLES = 600
LS_SMALL = [2, 4, 8, 16,32]
LS_LARGE = [32,64, 128, 256,512]
SEED = 123


def localizar_pasta_simulacao():
    pasta_atual = os.path.abspath(os.getcwd())
    if glob.glob(os.path.join(pasta_atual, "ciclo_*")):
        return pasta_atual
    nomes = {"simulação".casefold(), "simulaÃ§Ã£o".casefold()}

    if os.path.basename(pasta_atual).casefold() in nomes:
        return pasta_atual

    for nome in ("Simulação", "SimulaÃ§Ã£o"):
        pasta_direta = os.path.join(pasta_atual, nome)
        if os.path.isdir(pasta_direta):
            return pasta_direta

    candidatas = sorted(
        caminho
        for caminho in glob.glob(
            os.path.join(pasta_atual, "**", "Simula*"), recursive=True
        )
        if os.path.isdir(caminho)
        and os.path.basename(caminho).casefold().startswith("simula")
    )
    if not candidatas:
        raise RuntimeError(
            f"Nenhuma pasta 'Simulação' encontrada dentro de: {pasta_atual}"
        )
    if len(candidatas) > 1:
        print(
            f"Aviso: foram encontradas {len(candidatas)} pastas Simulação; "
            f"usando a primeira: {candidatas[0]}"
        )
    return candidatas[0]


BASE_DIR = localizar_pasta_simulacao()
PASTA_SAIDA = os.path.join(BASE_DIR, "SAR")
os.makedirs(PASTA_SAIDA, exist_ok=True)
print(f"Pasta Simulação encontrada: {BASE_DIR}")


def extrair_valor_da_pasta(caminho, prefixo):
    return os.path.basename(caminho).replace(prefixo, "").replace("=", "")


def ordenar_numericamente(lista):
    return sorted(lista, key=float)


def encontrar_configuracoes(pastas_ciclo):
    por_wh = {}
    for ciclo in pastas_ciclo:
        for arquivo in glob.glob(os.path.join(ciclo, "W*", "h*", "NRK=*", ARQUIVO_ALVO)):
            nrk = os.path.dirname(arquivo)
            h = os.path.dirname(nrk)
            w = os.path.dirname(h)
            chave = (extrair_valor_da_pasta(w, "W"), extrair_valor_da_pasta(h, "h"))
            valor_nrk = extrair_valor_da_pasta(nrk, "NRK=")
            if float(valor_nrk) not in (75, 100):
                por_wh.setdefault(chave, set()).add(valor_nrk)
    return [(w, h, ordenar_numericamente(nrks)) for (w, h), nrks in sorted(
        por_wh.items(), key=lambda item: (float(item[0][0]), float(item[0][1]))
    )]


def lista_para_grade(caminho, L_fixo=512):
    data = np.loadtxt(caminho, dtype=np.int64)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    idx = data[:, 0]
    sp = data[:, 1]

    N = int(idx.max())
    L = int(round(math.sqrt(N)))
    if L != L_fixo:
        L = L_fixo
        N = L * L

    grid_flat = np.zeros(N, dtype=np.int64)
    ok = (idx >= 1) & (idx <= N)
    grid_flat[idx[ok] - 1] = sp[ok]
    return grid_flat.reshape(L, L)


def riqueza_janela_periodica(grid, top, left, l, ignorar_zero=True):
    L = grid.shape[0]
    rows = np.arange(top, top + l) % L
    cols = np.arange(left, left + l) % L
    especies = np.unique(grid[np.ix_(rows, cols)])
    if ignorar_zero:
        especies = especies[especies != 0]
    return especies.size


def S_de_A_periodico(grid, l, n_samples=200, ignorar_zero=True, rng=None):
    L = grid.shape[0]
    if rng is None:
        rng = np.random.default_rng()

    tops = rng.integers(0, L, size=n_samples)
    lefts = rng.integers(0, L, size=n_samples)
    valores = np.empty(n_samples, dtype=float)
    for k in range(n_samples):
        valores[k] = riqueza_janela_periodica(
            grid, int(tops[k]), int(lefts[k]), l, ignorar_zero
        )
    return float(valores.mean()), float(valores.std(ddof=1))


def ajusta_z(A, S, mask):
    A = np.asarray(A, dtype=float)[mask]
    S = np.asarray(S, dtype=float)[mask]
    ok = (A > 0) & (S > 0)
    if np.sum(ok) < 2:
        return np.nan, np.nan

    z, intercepto = np.polyfit(np.log10(A[ok]), np.log10(S[ok]), 1)
    return float(z), float(10**intercepto)


def media_e_erro(valores):
    valores = np.asarray(valores, dtype=float)
    media = float(np.nanmean(valores))
    erro = float(np.nanstd(valores, ddof=1)) if len(valores) > 1 else 0.0
    return media, erro


def calcular_curva(caminho):
    grid = lista_para_grade(caminho, L_fixo=L_FIXO)
    rng = np.random.default_rng(SEED)
    return np.array(
        [
            S_de_A_periodico(
                grid,
                l,
                n_samples=N_SAMPLES,
                ignorar_zero=IGNORAR_ZERO,
                rng=rng,
            )[0]
            for l in LS
        ],
        dtype=float,
    )


def gerar_figura(W, H, NRKs, pastas_ciclo):
    A = np.array([l * l for l in LS], dtype=float)
    mask_small = np.isin(LS, LS_SMALL)
    mask_large = np.isin(LS, LS_LARGE)

    fig = plt.figure(figsize=(16, 7.2), dpi=180)
    grade = fig.add_gridspec(
        2,
        2,
        width_ratios=[1.65, 1.0],
        height_ratios=[1, 1],
        wspace=0.24,
        hspace=0.18,
    )
    ax_sar = fig.add_subplot(grade[:, 0])
    ax_z_large = fig.add_subplot(grade[0, 1])
    ax_z_small = fig.add_subplot(grade[1, 1], sharex=ax_z_large)
    resultados = []
    for idx, NRK in enumerate(NRKs):
        curvas = []
        for ciclo in pastas_ciclo:
            caminho = os.path.join(
                ciclo, f"W{W}", f"h{H}", f"NRK={NRK}", ARQUIVO_ALVO
            )
            if not os.path.exists(caminho):
                print(f"    Não encontrado: {caminho}")
                continue
            try:
                curvas.append(calcular_curva(caminho))
            except (OSError, ValueError, IndexError) as erro:
                print(f"    Não foi possível processar {caminho}: {erro}")

        if len(curvas) < 2:
            print(f"    Nenhum ciclo válido para NRK={NRK}")
            continue

        curvas = np.vstack(curvas)
        curva_media = curvas.mean(axis=0)
        curva_erro = (
            curvas.std(axis=0, ddof=1)
            if curvas.shape[0] > 1
            else np.zeros_like(curva_media)
        )
        z_small = [ajusta_z(A, curva, mask_small)[0] for curva in curvas]
        z_large = [ajusta_z(A, curva, mask_large)[0] for curva in curvas]
        z_small_media, z_small_erro = media_e_erro(z_small)
        z_large_media, z_large_erro = media_e_erro(z_large)

        resultados.append(
            {
                "nrk": float(NRK),
                "curva": curva_media,
                "z_small": z_small_media,
                "z_small_erro": z_small_erro,
                "z_large": z_large_media,
                "z_large_erro": z_large_erro,
            }
        )
        ax_sar.errorbar(
            A,
            curva_media,
            yerr=curva_erro,
            fmt="o",
            linestyle="none",
            color=cor_n(NRK),
            markersize=6,
            elinewidth=1.1,
            capsize=2.5,
            alpha=0.9,
            label=rf"$n_r={NRK}$",
        )
        xs, ys = curva_suave(A, curva_media, log_x=True, log_y=True)
        ax_sar.plot(xs, ys, color=cor_n(NRK), linewidth=2.0)
        print(
            f"  NRK={NRK}: ciclos={len(curvas)} | "
            f"z_small={z_small_media:.4f} ± {z_small_erro:.4f} | "
            f"z_large={z_large_media:.4f} ± {z_large_erro:.4f}"
        )

    if not resultados:
        plt.close(fig)
        print(f"Nenhum resultado para W={W} | h={H}")
        return

    resultados.sort(key=lambda item: item["nrk"])
    nrks = np.array([item["nrk"] for item in resultados])
    curvas_medias = np.vstack([item["curva"] for item in resultados])

    # Os traços representam os dois regimes do ajuste piecewise, calculados
    # sobre a curva média entre todos os valores de n.
    curva_referencia = curvas_medias.mean(axis=0)
    z_ref_small, c_ref_small = ajusta_z(A, curva_referencia, mask_small)
    z_ref_large, c_ref_large = ajusta_z(A, curva_referencia, mask_large)

    A_small = A[mask_small]
    x_small = np.logspace(np.log10(A_small.min()), np.log10(A_small.max()), 100)
    ax_sar.plot(
        x_small,
        c_ref_small * x_small**z_ref_small,
        color="#555555",
        linestyle="--",
        linewidth=1.8,
        label="Ajuste — pequenas escalas",
    )
    A_large = A[mask_large]
    x_large = np.logspace(np.log10(A_large.min()), np.log10(A_large.max()), 100)
    ax_sar.plot(
        x_large,
        c_ref_large * x_large**z_ref_large,
        color="#888888",
        linestyle=":",
        linewidth=2.0,
        label="Ajuste — grandes escalas",
    )

    z_large = np.array([item["z_large"] for item in resultados])
    z_small = np.array([item["z_small"] for item in resultados])
    ax_z_large.errorbar(
        nrks,
        z_large,
        yerr=[item["z_large_erro"] for item in resultados],
        fmt="s",
        linestyle="none",
        color=cor_n(3),
        linewidth=1.4,
        markersize=6,
        capsize=4,
    )
    xz, yz = curva_suave(nrks, z_large)
    ax_z_large.plot(xz, yz, color=cor_n(3), linewidth=2.0)
    ax_z_small.errorbar(
        nrks,
        z_small,
        yerr=[item["z_small_erro"] for item in resultados],
        fmt="o",
        linestyle="none",
        color=cor_n(25),
        linewidth=1.4,
        markersize=6,
        capsize=4,
    )
    xz, yz = curva_suave(nrks, z_small)
    ax_z_small.plot(xz, yz, color=cor_n(25), linewidth=2.0)

    ax_sar.set_xscale("log")
    ax_sar.set_yscale("log")
    # Destaca visualmente os mesmos intervalos usados nos dois ajustes.
    ax_sar.axvspan(A_small.min(), A_small.max(), color="#dbe9f6", alpha=0.48,
                   zorder=-10)
    ax_sar.axvspan(A_large.min(), A_large.max(), color="#eef5fa", alpha=0.62,
                   zorder=-10)
    ax_sar.axvline(A_small.max(), color="#7795ad", linestyle="--",
                   linewidth=1.0, alpha=0.65, zorder=-5)
    ax_sar.axvline(A_large.min(), color="#7795ad", linestyle="--",
                   linewidth=1.0, alpha=0.65, zorder=-5)
    ax_sar.text(np.sqrt(A_small.min() * A_small.max()), 0.035,
                "Pequenas escalas", transform=ax_sar.get_xaxis_transform(),
                ha="center", va="bottom", color="#355d7a", fontsize=10)
    ax_sar.text(np.sqrt(A_large.min() * A_large.max()), 0.035,
                "Grandes escalas", transform=ax_sar.get_xaxis_transform(),
                ha="center", va="bottom", color="#5e7d94", fontsize=10)
    ax_sar.text(np.sqrt(A_small.max() * A_large.min()), 0.035,
                "transição", transform=ax_sar.get_xaxis_transform(),
                ha="center", va="bottom", color="#777777", fontsize=8,
                rotation=90)
    ax_sar.set_xlabel("A")
    ax_sar.set_ylabel(r"$\langle S \rangle$")
    ax_sar.grid(True, which="major", linestyle="--", linewidth=0.7, alpha=0.28)
    handles, labels = ax_sar.get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(.5,.01),
               frameon=False, ncol=max(1, len(labels)), fontsize=9,
               columnspacing=1.2, handlelength=2.4)

    ax_z_large.set_ylabel("z")
    ax_z_large.text(0.03, 0.92, "Grandes escalas", transform=ax_z_large.transAxes,
                    va="top", color=cor_n(3), fontsize=11)
    ax_z_large.tick_params(labelbottom=False)
    ax_z_small.set_xlabel(r"$n_r$")
    ax_z_small.set_ylabel("z")
    ax_z_small.text(0.03, 0.92, "Pequenas escalas", transform=ax_z_small.transAxes,
                    va="top", color=cor_n(25), fontsize=11)
    for eixo in (ax_z_large, ax_z_small):
        eixo.set_xticks(nrks)
        eixo.grid(True, linestyle="--", linewidth=0.7, alpha=0.28)
        eixo.margins(x=0.08, y=0.16)

    nome = f"W{W}_h{H}_SAR_MODELO.png"
    caminho_saida = os.path.join(PASTA_SAIDA, nome)
    fig.subplots_adjust(bottom=.13)
    plt.savefig(caminho_saida, bbox_inches="tight", dpi=170)
    plt.close(fig)
    print(f"Salvo: {caminho_saida}")


def gerar_curvas_sar_por_h(W,Hs,NRKs,dados,mapa_por_h,mp,mg,A):
    """Salva uma figura de curvas SAR para cada valor de H."""
    from matplotlib.lines import Line2D

    marcadores=("o","s","^","D","v","P","X","<",">")
    for H in Hs:
        nrks_h=[n for n in NRKs if (H,n) in dados]
        if not nrks_h:
            print(f"W={W} H={H}: nenhuma configuração válida para as curvas SAR")
            continue

        niveis=np.linspace(.42,.92,len(nrks_h))
        fig,ax=plt.subplots(figsize=(9.7,5.9),dpi=180)
        legenda_n=[]

        for i,n in enumerate(nrks_h):
            media,erro,_,_=dados[(H,n)]
            cor=mapa_por_h[H](niveis[i])
            marcador=marcadores[i%len(marcadores)]
            ax.errorbar(A,media,yerr=erro,fmt=marcador,ls="none",color=cor,
                        ms=4.0,elinewidth=.7,capsize=1.5,alpha=.75,zorder=6)
            for mask,linha in ((mp,"--"),(mg,"-")):
                z,c=ajusta_z(A,media,mask)
                xf=np.logspace(np.log10(A[mask].min()),np.log10(A[mask].max()),120)
                ax.plot(xf,c*xf**z,color=cor,ls=linha,lw=1.65,alpha=.82,zorder=5)
            legenda_n.append(Line2D([0],[0],color=cor,marker=marcador,lw=1.8,
                                    label=rf"$n_r={n}$"))

        ax.set_xscale("log"); ax.set_yscale("log")
        ax.grid(True,which="major",ls="--",lw=.55,alpha=.25)
        ax.set_xlabel(r"$A$"); ax.set_ylabel(r"$\langle S\rangle$")
        fig.legend(handles=legenda_n,loc="lower center",bbox_to_anchor=(.5,.075),
                   ncol=len(legenda_n),frameon=False,columnspacing=1.2,
                   handletextpad=.45)
        legenda_escala=[
            Line2D([0],[0],color="#444444",ls="--",lw=2,
                   label="Small scales"),
            Line2D([0],[0],color="#444444",ls="-",lw=2,
                   label="Large scales"),
        ]
        fig.legend(handles=legenda_escala,loc="lower center",ncol=2,
                   bbox_to_anchor=(.5,.012),frameon=False,handlelength=2.7,
                   columnspacing=2.6)
        fig.subplots_adjust(left=.10,right=.98,top=.97,bottom=.25)
        arq=os.path.join(PASTA_SAIDA,f"W{W}_h{H}_SAR_CURVAS.png")
        fig.savefig(arq,bbox_inches="tight",dpi=180)
        plt.close(fig)
        print(f"Salvo: {arq}")


def gerar_paineis(W, Hs, NRKs, ciclos):
    A=np.array([l*l for l in LS],float); mp=np.isin(LS,LS_SMALL); mg=np.isin(LS,LS_LARGE)
    dados={}
    for H in Hs:
        for n in NRKs:
            curvas=[]
            for ciclo in ciclos:
                arq=os.path.join(ciclo,f"W{W}",f"h{H}",f"NRK={n}",ARQUIVO_ALVO)
                if not os.path.exists(arq): continue
                try: curvas.append(calcular_curva(arq))
                except (OSError,ValueError,IndexError) as erro: print(f"Falha em {arq}: {erro}")
            if len(curvas)<2: continue
            curvas=np.vstack(curvas)
            zp=media_e_erro([ajusta_z(A,c,mp)[0] for c in curvas])
            zg=media_e_erro([ajusta_z(A,c,mg)[0] for c in curvas])
            dados[(H,n)]=(curvas.mean(0),curvas.std(0,ddof=1),zp,zg)
            print(f"W={W} h={H} n={n}: {len(curvas)} ciclos")

    NRKs=[n for n in NRKs if any((H,n) in dados for H in Hs)]
    if not NRKs:
        print(f"W={W}: nenhuma configuração válida para o SAR")
        return
    mapas_h=(plt.cm.Blues,plt.cm.Oranges,plt.cm.Greens,plt.cm.Purples,plt.cm.Reds)
    mapa_por_h={H:mapas_h[i%len(mapas_h)] for i,H in enumerate(Hs)}
    niveis_n=np.linspace(.42,.92,len(NRKs))
    marcadores=("o","s","^","D","v","P","X","<",">")
    indice_n={n:i for i,n in enumerate(NRKs)}
    fig,ax=plt.subplots(figsize=(9.7,5.9),dpi=180)
    encontrou=False
    for H in Hs:
        for n in NRKs:
            item=dados.get((H,n))
            if item is None: continue
            encontrou=True
            media,erro,_,_=item
            i=indice_n[n]
            cor=mapa_por_h[H](niveis_n[i])
            marcador=marcadores[i%len(marcadores)]
            ax.errorbar(A,media,yerr=erro,fmt=marcador,ls="none",color=cor,
                        ms=4.0,elinewidth=.7,capsize=1.5,alpha=.75,zorder=6)
            for mask,linha in ((mp,"--"),(mg,"-")):
                z,c=ajusta_z(A,media,mask)
                xf=np.logspace(np.log10(A[mask].min()),np.log10(A[mask].max()),120)
                ax.plot(xf,c*xf**z,color=cor,ls=linha,lw=1.65,alpha=.78,zorder=5)
    if not encontrou: ax.text(.5,.5,"sem dados",transform=ax.transAxes,ha="center")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.grid(True,which="major",ls="--",lw=.55,alpha=.25)
    ax.set_xlabel(r"$A$"); ax.set_ylabel(r"$\langle S\rangle$")
    from matplotlib.lines import Line2D
    legenda_h_curvas=[Line2D([0],[0],color=mapa_por_h[H](.72),lw=2.6,
                             label=rf"$H={H}$") for H in Hs]
    fig.legend(handles=legenda_h_curvas,loc="lower center",ncol=len(Hs),
               bbox_to_anchor=(.5,.135),frameon=False,columnspacing=2.0)
    legenda_n=[Line2D([0],[0],color=plt.cm.Greys(niveis_n[indice_n[n]]),
                      marker=marcadores[indice_n[n]%len(marcadores)],lw=1.8,
                      label=rf"$n_r={n}$") for n in NRKs]
    fig.legend(handles=legenda_n,loc="lower center",bbox_to_anchor=(.5,.075),
               ncol=len(NRKs),frameon=False,columnspacing=1.2,handletextpad=.45)
    legenda_escala=[Line2D([0],[0],color="#444444",ls="--",lw=2,label="Small scales"),
                    Line2D([0],[0],color="#444444",ls="-",lw=2,label="Large scales")]
    fig.legend(handles=legenda_escala,loc="lower center",ncol=2,
               bbox_to_anchor=(.5,.012),frameon=False,handlelength=2.7,
               columnspacing=2.6)
    fig.subplots_adjust(left=.10,right=.98,top=.97,bottom=.30)
    # O painel combinado deixou de ser salvo. As curvas são exportadas abaixo,
    # em uma figura independente para cada valor de H.
    plt.close(fig)
    gerar_curvas_sar_por_h(W,Hs,NRKs,dados,mapa_por_h,mp,mg,A)
    painel_antigo=os.path.join(PASTA_SAIDA,f"W{W}_SAR_CURVAS_PAINEL.png")
    if os.path.exists(painel_antigo):
        os.remove(painel_antigo)
        print(f"Removido resultado combinado antigo: {painel_antigo}")

    fig,ax=plt.subplots(figsize=(7.4,4.8),dpi=180)
    paleta_h=("#2166ac","#d6604d","#1b9e77","#984ea3","#e6ab02")
    cores_h={H:paleta_h[i%len(paleta_h)] for i,H in enumerate(Hs)}
    estilos=((2,r"$z_s(n_r)$","-","o"),(3,r"$z_l(n_r)$","--","s"))
    for H in Hs:
        cor=cores_h[H]
        for indice,_,linha,marcador in estilos:
            pontos=[(float(n),*dados[(H,n)][indice]) for n in NRKs if (H,n) in dados]
            x=np.array([p[0] for p in pontos]); y=np.array([p[1] for p in pontos]); e=np.array([p[2] for p in pontos])
            xs,ys=curva_suave(x,y)
            ax.plot(xs,ys,color=cor,ls=linha,lw=2.1)
            ax.errorbar(x,y,yerr=e,fmt=marcador,ls="none",color=cor,
                        ms=5,elinewidth=1.1,capsize=3)
    ax.set_xticks([float(n) for n in NRKs]); ax.grid(True,ls="--",lw=.6,alpha=.28)
    ax.margins(x=.06,y=.13); ax.set_xlabel(r"$n_r$"); ax.set_ylabel(r"$z$")
    legenda_h=[Line2D([0],[0],color=cores_h[H],lw=2.5,label=rf"$H={H}$") for H in Hs]
    fig.legend(handles=legenda_h,loc="lower center",ncol=len(Hs),frameon=False,
               bbox_to_anchor=(.5,.075),columnspacing=2.0)
    legenda_z=[Line2D([0],[0],color="#333333",ls=linha,marker=marcador,
                      lw=2,label=rotulo) for _,rotulo,linha,marcador in estilos]
    fig.legend(handles=legenda_z,loc="lower center",ncol=2,frameon=False,
               bbox_to_anchor=(.5,.008),columnspacing=2.8)
    fig.subplots_adjust(left=.11,right=.97,top=.96,bottom=.25)
    arq=os.path.join(PASTA_SAIDA,f"W{W}_SAR_EXPOENTES_PAINEL.png")
    fig.savefig(arq,bbox_inches="tight",dpi=180); plt.close(fig); print(f"Salvo: {arq}")


pastas_ciclo = sorted(glob.glob(os.path.join(BASE_DIR, "ciclo_*")))
print(f"{len(pastas_ciclo)} ciclos encontrados")
if not pastas_ciclo:
    raise RuntimeError(f"Nenhuma pasta ciclo_* encontrada dentro de: {BASE_DIR}")

configs=encontrar_configuracoes(pastas_ciclo)
for W in ordenar_numericamente({w for w,_,_ in configs}):
    Hs=ordenar_numericamente({h for w,h,_ in configs if w==W})
    NRKs=ordenar_numericamente({n for w,_,ns in configs if w==W for n in ns})
    gerar_paineis(W,Hs,NRKs,pastas_ciclo)

for W, H, NRKs in []:
    print(f"\nProcessando SAR — W={W} | h={H}")
    gerar_figura(W, H, NRKs, pastas_ciclo)

print("\nSAR finalizada com sucesso!")
