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
  py lambda_py_fast.py --test 40
  py lambda_py_fast.py --pubkey 02... --lo 0 --bits 50 [--workers 8] [--batch 256]

AVISO HONESTO: o custo do Lambda e ~2*sqrt(2^bits) passos. Para 160 bits isso
e ~2^81 passos - inviavel em qualquer hardware existente (Python, C++ ou GPU).
O codigo aceita 160 bits, mas na pratica so e viavel ate ~55-60 bits em CPU.
"""
import argparse
import math
import multiprocessing as mp
import os
import random
import sys
import time
from queue import Empty

try:
    import gmpy2
    from gmpy2 import mpz, invert as _gmp_invert
    HAVE_GMPY = True
except ImportError:  # sem gmpy2: usa int nativo
    HAVE_GMPY = False

P = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F
N = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
G = (0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798,
     0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8)


# ----------------------------------------------------------------- curva
def add(a, b):
    """Soma de pontos (afim). None = infinito."""
    if a is None:
        return b
    if b is None:
        return a
    x1, y1 = a
    x2, y2 = b
    if x1 == x2:
        if (y1 + y2) % P == 0:
            return None
        lam = 3 * x1 * x1 * pow(2 * y1, -1, P) % P
    else:
        lam = (y2 - y1) * pow(x2 - x1, -1, P) % P
    x3 = (lam * lam - x1 - x2) % P
    return (x3, (lam * (x1 - x3) - y1) % P)


def mul(k, pt=G):
    k %= N
    r = None
    while k:
        if k & 1:
            r = add(r, pt)
        pt = add(pt, pt)
        k >>= 1
    return r


def parse_pubkey(h):
    h = h.strip().lower()
    if h.startswith("04") and len(h) == 130:
        return (int(h[2:66], 16), int(h[66:], 16))
    x = int(h[2:], 16)
    y = pow((x ** 3 + 7) % P, (P + 1) // 4, P)
    if (y & 1) != (int(h[:2], 16) & 1):
        y = P - y
    return (x, y)


# ------------------------------------------------------------ parametros
def pick_params(bits, total_kangaroos):
    """Tabela de saltos (potencias de 2) e tamanho do ponto distinto."""
    W = 1 << bits
    target = math.sqrt(W) / 2
    m = min(range(2, bits + 3), key=lambda m: abs((2 ** m - 1) / m - target))
    # desperdicio ~ kangaroos * 2^dbits deve ser << sqrt(W)
    dbits = max(0, bits // 2 - 4 - int(math.log2(max(1, total_kangaroos))))
    return m, dbits


# ---------------------------------------------------------------- worker
def worker(Q, lo, bits, K, m, dbits, queue, stop):
    random.seed(os.urandom(16))
    W = 1 << bits
    mid = lo + W // 2
    dmask = (1 << dbits) - 1
    limit = max(2000, 40 << dbits)  # passos sem DP => provavel ciclo, recomeca

    conv = mpz if HAVE_GMPY else int
    Pm = conv(P)
    if HAVE_GMPY:
        inv = lambda a: _gmp_invert(a, Pm)
    else:
        inv = lambda a: pow(a, -1, P)

    jumps = [1 << i for i in range(m)]
    jpts = [mul(s) for s in jumps]
    jx = [conv(p[0]) for p in jpts]
    jy = [conv(p[1]) for p in jpts]

    def start(kind):
        while True:
            if kind == 0:
                d = mid + random.randrange(-W // 8, W // 8 + 1)
                pt = mul(d)
            else:
                d = random.randrange(0, W // 8 + 1)
                pt = add(Q, mul(d))
            if pt is not None:
                return conv(pt[0]), conv(pt[1]), d

    X = [0] * K
    Y = [0] * K
    D = [0] * K
    kind = [i & 1 for i in range(K)]       # metade domados, metade selvagens
    cnt = [0] * K
    for i in range(K):
        X[i], Y[i], D[i] = start(kind[i])

    pref = [0] * K
    dxs = [0] * K
    steps_pending = 0
    batches = 0

    while not stop.is_set():
        idx = [int(x % m) for x in X]
        bad = []
        acc = conv(1)
        for i in range(K):
            d = (jx[idx[i]] - X[i]) % Pm
            if d == 0:
                d = conv(1)
                bad.append(i)
            dxs[i] = d
            pref[i] = acc
            acc = acc * d % Pm
        iv = inv(acc)
        out = []
        for i in range(K - 1, -1, -1):
            j = idx[i]
            d = dxs[i]
            invd = iv * pref[i] % Pm
            iv = iv * d % Pm
            x = X[i]
            lam = (jy[j] - Y[i]) * invd % Pm
            x3 = (lam * lam - x - jx[j]) % Pm
            Y[i] = (lam * (x - x3) - Y[i]) % Pm
            X[i] = x3
            D[i] += jumps[j]
            if (x3 & dmask) == 0:
                out.append((int(x3), D[i], kind[i]))
                cnt[i] = 0
            else:
                cnt[i] += 1
                if cnt[i] > limit:
                    bad.append(i)
        for i in bad:
            X[i], Y[i], D[i] = start(kind[i])
            cnt[i] = 0

        steps_pending += K
        batches += 1
        if out or batches % 20 == 0:
            queue.put((out, steps_pending))
            steps_pending = 0


# ----------------------------------------------------------- coordenador
def kangaroo(Q, lo, bits, workers=None, K=256, log=True):
    """Procura k em [lo, lo+2^bits) com G*k = Q."""
    workers = workers or os.cpu_count() or 1
    W = 1 << bits
    m, dbits = pick_params(bits, K * workers)
    if log:
        exp = 2.1 * math.sqrt(W)
        print(f"  workers={workers} canguros/worker={K} DP bits={dbits} "
              f"tabela de saltos={m} gmpy2={'sim' if HAVE_GMPY else 'nao'}")
        print(f"  passos esperados ~ 2^{math.log2(exp):.1f}")
        if bits > 60:
            print("  AVISO: acima de ~60 bits isso nao termina em tempo util.")

    queue = mp.Queue()
    stop = mp.Event()
    procs = [mp.Process(target=worker, args=(Q, lo, bits, K, m, dbits, queue, stop),
                        daemon=True) for _ in range(workers)]
    for p in procs:
        p.start()

    tables = ({}, {})  # 0 = domado, 1 = selvagem
    total = 0
    dps = 0
    t0 = time.time()
    last = t0
    result = None
    try:
        while result is None:
            try:
                out, steps = queue.get(timeout=1.0)
            except Empty:
                out, steps = [], 0
            total += steps
            for x, dist, kd in out:
                dps += 1
                other = tables[1 - kd].get(x)
                if other is not None:
                    k = (dist - other) if kd == 0 else (other - dist)
                    if 0 <= k < N and mul(k) == Q:
                        result = k
                        break
                tables[kd][x] = dist
            now = time.time()
            if log and now - last >= 5:
                last = now
                print(f"  {total:,} passos  {total/(now-t0):,.0f} passos/s  "
                      f"DPs={dps:,}  {now-t0:.0f}s")
    finally:
        stop.set()
        for p in procs:
            p.terminate()
    if log:
        el = time.time() - t0
        print(f"  passos: {total:,}  tempo: {el:.1f}s  ({total/max(el,1e-9):,.0f} passos/s)")
    return result


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pubkey", help="chave publica (hex, comprimida ou nao)")
    ap.add_argument("--lo", default="0", help="inicio do intervalo (hex)")
    ap.add_argument("--bits", type=int, help="tamanho do intervalo em bits (ate 160)")
    ap.add_argument("--test", type=int, help="auto-teste com chave aleatoria de N bits")
    ap.add_argument("--workers", type=int, default=None, help="processos (padrao: nucleos)")
    ap.add_argument("--batch", type=int, default=256, help="canguros por processo")
    a = ap.parse_args()

    if a.test:
        bits = a.test
        lo = 1 << (bits - 1)
        k = random.randrange(lo, lo + (1 << (bits - 1)))
        Q = mul(k)
        print(f"Teste: chave secreta de {bits} bits, procurando em [2^{bits-1}, 2^{bits})")
        t = time.time()
        r = kangaroo(Q, lo, bits - 1, a.workers, a.batch)
        print(f"Achou: {hex(r)}  (correto: {r == k})  em {time.time()-t:.1f}s")
        return

    if not (a.pubkey and a.bits):
        ap.error("use --test N  ou  --pubkey HEX --lo HEX --bits N")
    if a.bits > 160:
        ap.error("--bits maximo: 160")
    Q = parse_pubkey(a.pubkey)
    lo = int(a.lo, 16)
    r = kangaroo(Q, lo, a.bits, a.workers, a.batch)
    print(f"Chave privada: {hex(r)}")


if __name__ == "__main__":
    mp.freeze_support()
    main()
