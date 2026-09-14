import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from estilo_graficos import cor_n, curva_suave

# =========================
# CONFIGURAÇÕES
# =========================
PASTA_ATUAL = os.getcwd()

if glob.glob(os.path.join(PASTA_ATUAL, "ciclo_*")):
    BASE_DIR = PASTA_ATUAL
elif os.path.basename(PASTA_ATUAL).casefold() == "simulação".casefold():
    BASE_DIR = PASTA_ATUAL
else:
    pasta_direta = os.path.join(PASTA_ATUAL, "Simulação")
    if os.path.isdir(pasta_direta):
        BASE_DIR = pasta_direta
    else:
        candidatas = sorted(
            caminho
            for caminho in glob.glob(
                os.path.join(PASTA_ATUAL, "**", "Simulação"), recursive=True
            )
            if os.path.isdir(caminho)
        )
        if not candidatas:
            raise RuntimeError(
                f"Nenhuma pasta 'Simulação' encontrada dentro de: {PASTA_ATUAL}"
            )
        BASE_DIR = candidatas[0]
        if len(candidatas) > 1:
            print(
                f"⚠️ Foram encontradas {len(candidatas)} pastas Simulação; "
                f"usando a primeira: {BASE_DIR}"
            )

print(f"📁 Pasta Simulação encontrada: {BASE_DIR}")
ARQUIVO_ALVO = "N_ESP.dat"

PASTA_MEDIAS = os.path.join(BASE_DIR, "MEDIAS")
os.makedirs(PASTA_MEDIAS, exist_ok=True)
PASTA_INDIVIDUAIS = os.path.join(BASE_DIR, "INDIVIDUAIS")
os.makedirs(PASTA_INDIVIDUAIS, exist_ok=True)

# =========================
# FUNÇÕES
# =========================
def ler_n_esp(caminho):
    dados = np.loadtxt(caminho)

    if dados.ndim == 1:
        t = np.array([dados[0]])
        n = np.array([dados[1]])
    else:
        t = dados[:, 0]
        n = dados[:, 1]

    return t, n


def gerar_curva_individual(caminho):
    try:
        t, n = ler_n_esp(caminho)
    except (OSError, ValueError, IndexError) as erro:
        print(f"⚠️ Não foi possível ler {caminho}: {erro}")
        return False

    if len(t) == 0:
        print(f"⚠️ Arquivo vazio: {caminho}")
        return False

    pasta = os.path.dirname(caminho)
    configuracao = os.path.relpath(pasta, BASE_DIR)

    plt.figure(figsize=(10, 7), dpi=150)
    ts, ns = curva_suave(t, n)
    plt.plot(ts, ns, color=cor_n(os.path.basename(pasta).replace("NRK=", "")), linewidth=1.8)
    plt.plot(t, n, "o", color=cor_n(os.path.basename(pasta).replace("NRK=", "")), markersize=2.5)
    plt.xlabel("Time")
    plt.ylabel("S")
    plt.grid(False)
    plt.tight_layout()

    pasta_saida = os.path.join(PASTA_INDIVIDUAIS, configuracao)
    os.makedirs(pasta_saida, exist_ok=True)
    saida = os.path.join(pasta_saida, "Esp_t.png")
    plt.savefig(saida, bbox_inches="tight")
    plt.close()
    print(f"✅ Curva individual: {saida}")
    return True


def extrair_valor_da_pasta(caminho, prefixo):
    base = os.path.basename(caminho)
    return base.replace(prefixo, "").replace("=", "")


def ordenar_numericamente(lista):
    return sorted(lista, key=lambda x: float(x))


def encontrar_configuracoes(pastas_ciclo):
    por_wh = {}
    for ciclo in pastas_ciclo:
        padrao = os.path.join(ciclo, "W*", "h*", "NRK=*", ARQUIVO_ALVO)
        for arquivo in glob.glob(padrao):
            pasta_nrk = os.path.dirname(arquivo)
            pasta_h = os.path.dirname(pasta_nrk)
            pasta_w = os.path.dirname(pasta_h)
            W = extrair_valor_da_pasta(pasta_w, "W")
            H = extrair_valor_da_pasta(pasta_h, "h")
            NRK = extrair_valor_da_pasta(pasta_nrk, "NRK=")
            if float(NRK) not in (75, 100):
                por_wh.setdefault((W, H), set()).add(NRK)
    return [
        (W, H, ordenar_numericamente(NRKs))
        for (W, H), NRKs in sorted(
            por_wh.items(), key=lambda item: (float(item[0][0]), float(item[0][1]))
        )
    ]


# =========================
# LOCALIZA CICLOS
# =========================
pastas_ciclo = sorted(glob.glob(os.path.join(BASE_DIR, "ciclo_*")))

print(f"\n🔍 {len(pastas_ciclo)} ciclos encontrados\n")

if len(pastas_ciclo) == 0:
    raise RuntimeError(f"Nenhuma pasta ciclo_* encontrada dentro de: {BASE_DIR}")

# Gera uma curva própria ao lado de cada N_ESP.dat.
arquivos_n_esp = sorted(
    glob.glob(os.path.join(BASE_DIR, "ciclo_*", "**", ARQUIVO_ALVO), recursive=True)
)
arquivos_n_esp = [
    arquivo
    for arquivo in arquivos_n_esp
    if os.path.basename(os.path.dirname(arquivo)) not in {"NRK=75", "NRK=100"}
]
print(f"📈 Gerando {len(arquivos_n_esp)} curvas individuais...")
for arquivo_n_esp in arquivos_n_esp:
    gerar_curva_individual(arquivo_n_esp)

configuracoes = encontrar_configuracoes(pastas_ciclo)

# =========================
# LOOP PRINCIPAL
# =========================
for W, H, NRKs in configuracoes:

    print(f"\n📊 Processando médias para W={W} | h={H}")

    plt.figure(figsize=(10, 7), dpi=150)

    for idx, NRK in enumerate(NRKs):

        series = []
        t_ref = None

        for ciclo in pastas_ciclo:

            caminho = os.path.join(
                ciclo,
                f"W{W}",
                f"h{H}",
                f"NRK={NRK}",
                ARQUIVO_ALVO
            )

            if not os.path.exists(caminho):
                print(f"⚠️ Não encontrado: {caminho}")
                continue

            t, n = ler_n_esp(caminho)

            if t_ref is None:
                t_ref = t

            series.append(n)

        if len(series) < 2:
            print(f"⚠️ Nenhuma série encontrada para NRK={NRK}")
            continue

        min_len = min(len(s) for s in series)
        series = np.array([s[:min_len] for s in series])

        t_plot = t_ref[:min_len]

        media = np.mean(series, axis=0)
        std = np.std(series, axis=0)

        ts, medias = curva_suave(t_plot, media)
        _, desvios = curva_suave(t_plot, std)
        plt.plot(
            ts,
            medias,
            color=cor_n(NRK),
            linewidth=1.8,
            label=rf"$n_r={NRK}$"
        )
        plt.plot(t_plot, media, "o", color=cor_n(NRK), markersize=2.2)

        plt.fill_between(
            ts,
            medias - desvios,
            medias + desvios,
            color=cor_n(NRK),
            alpha=0.25
        )

    plt.xlabel("Time")
    plt.ylabel(r"$\langle S \rangle$")

    plt.grid(False)
    plt.legend(loc="upper center", bbox_to_anchor=(0.5, -0.14),
               ncol=max(1, len(NRKs)), frameon=False)
    plt.subplots_adjust(bottom=0.22)

    nome = f"W{W}_h{H}_MEDIA.png"
    plt.savefig(os.path.join(PASTA_MEDIAS, nome), bbox_inches="tight")
    plt.close()

    print(f"✅ Salvo: {os.path.join('MEDIAS', nome)}")

print("\n🏁 MÉDIAS FINALIZADAS COM SUCESSO!")
