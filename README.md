#!/usr/bin/env python3
"""
Pollard's Lambda (kangaroo) para secp256k1 - versao rapida, Python 3.

Otimizacoes em relacao a lambda_py.py:
  1. Inversao em lote (truque de Montgomery): 1 inversao modular para cada
     K saltos, em vez de 1 por salto (a inversao e o gargalo no Python).
  2. Varios canguros por processo (K) avancando juntos em lote.
  3. Multiprocessamento: um processo por nucleo; os pontos distintos (DP)
     vao para um processo coordenador que detecta as colisoes.
  4. Usa gmpy2 automaticamente se estiver instalado (pip install gmpy2).
  5. Faixa de intervalo ate 160 bits aceita (inteiros arbitrarios).

Uso:
  py lambda_py_fast.py --pubkey 02... --lo 0 --bits 50 [--workers 8] [--batch 256]

CONFIGURAÇÃO

1. Instalar python 3.12
2. Baixar ZIP
3. cd
4.  py lambda_py_fast.py --pubkey 03a2efa402fd5268400c77c20e574ba86409ededee7c4020e4b9f0edbee53de0d4 --lo 8000000000 bits 40 --workers 8 --batch 256
5.  py lambda_py_fast.py --pubkey 031f6a332d3c5c4f2de2378c012f429cd109ba07d69690c6c701b6bb87860d6640 --lo 8000000000000000000000000000000000000 bits 140 --workers 8 --batch 256
