#!/usr/bin/env python3
"""Gera o numero de especies no ultimo passo em funcao de NRK.

Cada ponto individual usa a segunda coluna do ultimo registro de N_ESP.dat.
Graficos individuais aceitam dados incompletos; graficos medios exigem ao menos
dois ciclos para cada combinacao W/h/NRK.
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path
from statistics import fmean, stdev

import matplotlib.pyplot as plt
from estilo_graficos import cor_h, curva_suave


def numero(nome: str, prefixo: str) -> float:
    return float(nome.removeprefix(prefixo).removeprefix("="))


def especies_no_ultimo_passo(arquivo: Path) -> float:
    """Retorna S no registro com o maior passo de tempo de N_ESP.dat."""
    ultimo: tuple[float, float] | None = None
    with arquivo.open(encoding="utf-8", errors="ignore") as entrada:
        for linha_numero, linha in enumerate(entrada, start=1):
            colunas = linha.split()
            if not colunas:
                continue
            if len(colunas) < 2:
                raise ValueError(f"{arquivo}, linha {linha_numero}: esperadas 2 colunas")
            try:
                registro = (float(colunas[0]), float(colunas[1]))
            except ValueError as erro:
                raise ValueError(
                    f"{arquivo}, linha {linha_numero}: valor numerico invalido"
                ) from erro
            if ultimo is None or registro[0] >= ultimo[0]:
                ultimo = registro
    if ultimo is None:
        raise ValueError(f"Arquivo vazio: {arquivo}")
    return ultimo[1]


def coletar(base: Path):
    # chave: (W, h, NRK); valor: [(ciclo, especies_no_ultimo_passo)]
    dados = defaultdict(list)
    for ciclo in sorted(p for p in base.glob("ciclo_*") if p.is_dir()):
        for arquivo in ciclo.glob("W*/h*/NRK=*/N_ESP.dat"):
            nrk_dir = arquivo.parent
            h_dir = nrk_dir.parent
            w_dir = h_dir.parent
            chave = (
                numero(w_dir.name, "W"),
                numero(h_dir.name, "h"),
                numero(nrk_dir.name, "NRK"),
            )
            if chave[2] in (75, 100):
                continue
            try:
                valor = especies_no_ultimo_passo(arquivo)
            except (OSError, ValueError) as erro:
                print(f"AVISO: {erro}")
                continue
            dados[chave].append((ciclo.name, valor))
    return dict(dados)


def gerar_individuais(destino: Path, dados, dpi: int) -> int:
    por_ciclo_w = defaultdict(lambda: defaultdict(list))
    for (w, h, nrk), observacoes in dados.items():
        for ciclo, valor in observacoes:
            por_ciclo_w[(ciclo, w)][h].append((nrk, valor))

    gerados = 0
    for (ciclo, w), por_h in sorted(por_ciclo_w.items()):
        fig, ax = plt.subplots(figsize=(10, 6.5))
        valores_nrk = set()
        for h, pontos in sorted(por_h.items()):
            pontos.sort()
            valores_nrk.update(p[0] for p in pontos)
            x, y = [p[0] for p in pontos], [p[1] for p in pontos]
            xs, ys = curva_suave(x, y)
            ax.plot(xs, ys, color=cor_h(h), linewidth=1.8, label=f"H={h:.2f}")
            ax.plot(x, y, "o", color=cor_h(h), markersize=5)
        ax.set_xlabel(r"$n_r$")
        ax.set_ylabel(r"$\langle S \rangle$")
        ax.set_xticks(sorted(valores_nrk))
        ax.grid(True, linestyle="--", alpha=0.35)
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.14),
                  ncol=max(1, len(por_h)), frameon=False)
        fig.subplots_adjust(bottom=0.22)
        pasta = destino / "INDIVIDUAIS" / ciclo
        pasta.mkdir(parents=True, exist_ok=True)
        saida = pasta / f"populacao_x_NRK_W{w:g}.png"
        fig.savefig(saida, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"Individual: {saida}")
        gerados += 1
    return gerados


def gerar_medias(destino: Path, dados, dpi: int) -> tuple[int, list[list[object]]]:
    por_w = defaultdict(lambda: defaultdict(list))
    linhas_csv: list[list[object]] = []
    for (w, h, nrk), observacoes in sorted(dados.items()):
        if len(observacoes) < 2:
            print(f"Sem complemento: W={w:g}, h={h:.2f}, NRK={nrk:g}")
            continue
        valores = [valor for _, valor in observacoes]
        media = fmean(valores)
        desvio = stdev(valores)
        por_w[w][h].append((nrk, media, desvio))
        linhas_csv.append(
            [w, h, nrk, len(valores), media, desvio, ";".join(c for c, _ in observacoes)]
        )

    gerados = 0
    pasta = destino / "MEDIAS"
    pasta.mkdir(parents=True, exist_ok=True)
    for w, por_h in sorted(por_w.items()):
        fig, ax = plt.subplots(figsize=(10, 6.5))
        valores_nrk = set()
        for h, pontos in sorted(por_h.items()):
            pontos.sort()
            valores_nrk.update(p[0] for p in pontos)
            x, y = [p[0] for p in pontos], [p[1] for p in pontos]
            xs, ys = curva_suave(x, y)
            ax.plot(xs, ys, color=cor_h(h), linewidth=1.8, label=f"H={h:.2f}")
            ax.errorbar(
                x,
                y,
                yerr=[p[2] for p in pontos],
                marker="o",
                linestyle="none",
                color=cor_h(h),
                linewidth=1.2,
                elinewidth=1,
                capsize=3,
            )
        ax.set_xlabel(r"$n_r$")
        ax.set_ylabel(r"$\langle S \rangle$")
        ax.set_xticks(sorted(valores_nrk))
        ax.grid(True, linestyle="--", alpha=0.35)
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.14),
                  ncol=max(1, len(por_h)), frameon=False)
        fig.subplots_adjust(bottom=0.22)
        saida = pasta / f"populacao_x_NRK_W{w:g}_MEDIA.png"
        fig.savefig(saida, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"Media: {saida}")
        gerados += 1
    return gerados, linhas_csv


def salvar_csv(destino: Path, linhas: list[list[object]]) -> Path:
    saida = destino / "especies_ultimo_passo_x_NRK.csv"
    with saida.open("w", encoding="utf-8", newline="") as arquivo:
        escritor = csv.writer(arquivo)
        escritor.writerow(["W", "h", "NRK", "n_ciclos", "media", "desvio", "ciclos"])
        escritor.writerows(linhas)
    return saida


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--resultados",
        type=Path,
        default=Path(__file__).resolve().parents[2] / "dados" / "resultados",
    )
    parser.add_argument("--dpi", type=int, default=300)
    args = parser.parse_args()

    base = args.resultados.expanduser().resolve()
    if not base.is_dir():
        parser.error(f"pasta nao encontrada: {base}")
    destino = base / "POP_X_NRK"
    destino.mkdir(parents=True, exist_ok=True)

    dados = coletar(base)
    if not dados:
        print("Nenhum N_ESP.dat valido encontrado.")
        return 1
    individuais = gerar_individuais(destino, dados, args.dpi)
    medias, linhas = gerar_medias(destino, dados, args.dpi)
    tabela = salvar_csv(destino, linhas)
    print(f"Concluido: {individuais} individuais, {medias} medias; tabela: {tabela}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
