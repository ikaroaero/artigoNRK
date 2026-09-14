"""Executa todas as analises sobre as bases organizadas de resultados.

As analises medias usam apenas configuracoes presentes em pelo menos dois
ciclos. As curvas Esp_t individuais sao geradas para toda configuracao que
possua N_ESP.dat, mesmo sem ciclo complementar.

Sem argumentos, processa ``dados/resultados`` (W=4) e
``dados/importados_maquina58`` (W=3), mantendo as saidas dentro de cada base.
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

#ANALISES = ("SAR.py")

ANALISES = ("Esp_t.py", "SAR.py", "SAD.py", "SAD_preston.py")


def inventario(base: Path) -> dict[Path, list[Path]]:
    encontrados: dict[Path, list[Path]] = defaultdict(list)
    for ciclo in sorted(base.glob("ciclo_*")):
        for arquivo in ciclo.glob("W*/h*/NRK=*/N_ESP.dat"):
            chave = arquivo.parent.relative_to(ciclo)
            encontrados[chave].append(ciclo)
    return dict(encontrados)


def processar_base(base: Path, apenas_inventario: bool) -> int:
    """Inventaria e executa todas as analises para uma unica base."""
    if not base.is_dir():
        print(f"Pasta nao encontrada: {base}", file=sys.stderr)
        return 1

    dados = inventario(base)
    individuais = len(dados)
    complementares = {k: v for k, v in dados.items() if len(v) >= 2}
    sem_complemento = {k: v for k, v in dados.items() if len(v) < 2}

    print(f"Pasta: {base}")
    print(f"Configuracoes distintas: {individuais}")
    print(f"Com ciclo complementar: {len(complementares)}")
    print(f"Sem ciclo complementar: {len(sem_complemento)}")
    for chave, ciclos in sorted(dados.items(), key=lambda item: str(item[0])):
        nomes = ", ".join(ciclo.name for ciclo in ciclos)
        tipo = "MEDIA" if len(ciclos) >= 2 else "APENAS INDIVIDUAL"
        print(f"  [{tipo}] {chave} <- {nomes}")

    if apenas_inventario:
        return 0

    pasta_scripts = Path(__file__).resolve().parent
    ambiente = os.environ.copy()
    ambiente["PYTHONIOENCODING"] = "utf-8"
    for nome in ANALISES:
        script = pasta_scripts / nome
        print(f"\n===== Executando {nome} =====", flush=True)
        resultado = subprocess.run([sys.executable, str(script)], cwd=base, env=ambiente)
        if resultado.returncode != 0:
            print(f"Falha em {nome} (codigo {resultado.returncode}).", file=sys.stderr)
            return resultado.returncode

    print("\n===== Executando Pop_NRK.py =====", flush=True)
    resultado = subprocess.run(
        [sys.executable, str(pasta_scripts / "Pop_NRK.py"), "--resultados", str(base)],
        cwd=base,
        env=ambiente,
    )
    if resultado.returncode != 0:
        print(f"Falha em Pop_NRK.py (codigo {resultado.returncode}).", file=sys.stderr)
        return resultado.returncode

    print(f"\nAnalises concluidas para: {base}")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Gera Esp_t individual e analises medias entre ciclos complementares."
    )
    parser.add_argument(
        "--resultados",
        type=Path,
        action="append",
        help=(
            "Pasta que contem ciclo_001, ciclo_002 etc. Pode ser repetido. "
            "Se omitido, processa resultados (W=4) e importados_maquina58 (W=3)."
        ),
    )
    parser.add_argument(
        "--inventario",
        action="store_true",
        help="Apenas mostra os complementos encontrados, sem gerar analises.",
    )
    args = parser.parse_args()

    raiz = Path(__file__).resolve().parents[2]
    bases = args.resultados or [
        raiz / "dados" / "resultados",
        raiz / "dados" / "importados_maquina58",
    ]

    for indice, pasta in enumerate(bases):
        base = pasta.expanduser().resolve()
        if indice:
            print("\n" + "=" * 72)
        print(f"BASE DE ANALISE: {base}")
        codigo = processar_base(base, args.inventario)
        if codigo != 0:
            return codigo

    print("\nTodas as bases foram processadas com sucesso.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
